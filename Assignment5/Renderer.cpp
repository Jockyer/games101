#include <fstream>
#include "Vector.hpp"
#include "Renderer.hpp"
#include "Scene.hpp"
#include <optional>

inline float deg2rad(const float &deg)
{ return deg * M_PI/180.0; }

// Compute reflection direction
Vector3f reflect(const Vector3f &I, const Vector3f &N)
{
    return I - 2 * dotProduct(I, N) * N;
}

// [comment]
// Compute refraction direction using Snell's law
//
// We need to handle with care the two possible situations:
//
//    - When the ray is inside the object
//
//    - When the ray is outside.
//
// If the ray is outside, you need to make cosi positive cosi = -N.I
//
// If the ray is inside, you need to invert the refractive indices and negate the normal N
// [/comment]
Vector3f refract(const Vector3f &I, const Vector3f &N, const float &ior)
{
    float cosi = clamp(-1, 1, dotProduct(I, N));
    float etai = 1, etat = ior;
    Vector3f n = N;
    if (cosi < 0) { cosi = -cosi; } else { std::swap(etai, etat); n= -N; }
    float eta = etai / etat;
    float k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? 0 : eta * I + (eta * cosi - sqrtf(k)) * n;
}

// [comment]
// Compute Fresnel equation
//
// \param I is the incident view direction
//
// \param N is the normal at the intersection point
//
// \param ior is the material refractive index
// [/comment]
float fresnel(const Vector3f &I, const Vector3f &N, const float &ior)
{
    float cosi = clamp(-1, 1, dotProduct(I, N));
    float etai = 1, etat = ior;
    if (cosi > 0) {  std::swap(etai, etat); }
    // Compute sini using Snell's law
    float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
    // Total internal reflection
    if (sint >= 1) {
        return 1;
    }
    else {
        float cost = sqrtf(std::max(0.f, 1 - sint * sint));
        cosi = fabsf(cosi);
        float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) / 2;
    }
    // As a consequence of the conservation of energy, transmittance is given by:
    // kt = 1 - kr;
}

// [comment]
// Returns true if the ray intersects an object, false otherwise.
//
// \param orig is the ray origin
// \param dir is the ray direction
// \param objects is the list of objects the scene contains
// \param[out] tNear contains the distance to the cloesest intersected object.
// \param[out] index stores the index of the intersect triangle if the interesected object is a mesh.
// \param[out] uv stores the u and v barycentric coordinates of the intersected point
// \param[out] *hitObject stores the pointer to the intersected object (used to retrieve material information, etc.)
// \param isShadowRay is it a shadow ray. We can return from the function sooner as soon as we have found a hit.
// [/comment]
std::optional<hit_payload> trace(
        const Vector3f &orig, const Vector3f &dir,
        const std::vector<std::unique_ptr<Object> > &objects)
{
    float tNear = kInfinity;
    std::optional<hit_payload> payload;
    for (const auto & object : objects)
    {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (object->intersect(orig, dir, tNearK, indexK, uvK) && tNearK < tNear)
        {
            payload.emplace();
            payload->hit_obj = object.get();
            payload->tNear = tNearK;
            payload->index = indexK;
            payload->uv = uvK;
            tNear = tNearK;
        }
    }

    return payload;
}

// [comment]
// Implementation of the Whitted-style light transport algorithm (E [S*] (D|G) L)
//
// This function is the function that compute the color at the intersection point
// of a ray defined by a position and a direction. Note that thus function is recursive (it calls itself).
//
// If the material of the intersected object is either reflective or reflective and refractive,
// then we compute the reflection/refraction direction and cast two new rays into the scene
// by calling the castRay() function recursively. When the surface is transparent, we mix
// the reflection and refraction color using the result of the fresnel equations (it computes
// the amount of reflection and refraction depending on the surface normal, incident view direction
// and surface refractive index).
//
// If the surface is diffuse/glossy we use the Phong illumation model to compute the color
// at the intersection point.
// [/comment]
Vector3f castRay(
        const Vector3f &orig, const Vector3f &dir, const Scene& scene,
        int depth)
{
    if (depth > scene.maxDepth) {
        return Vector3f(0.0,0.0,0.0);
    }

    Vector3f hitColor = scene.backgroundColor;
    if (auto payload = trace(orig, dir, scene.get_objects()); payload)
    {
        Vector3f hitPoint = orig + dir * payload->tNear;
        Vector3f N; // normal
        Vector2f st; // st coordinates
        payload->hit_obj->getSurfaceProperties(hitPoint, dir, payload->index, payload->uv, N, st);
        switch (payload->hit_obj->materialType) {
            case REFLECTION_AND_REFRACTION:
            {
                // TODO: 反射/折射两条递归光线 + fresnel 混合
                float kr = fresnel(dir, N, payload->hit_obj->ior);
                Vector3f reflDir = reflect(dir, N);
                Vector3f refrDir = refract(dir, N, payload->hit_obj->ior);
                
                reflDir = normalize(reflDir);
                refrDir = normalize(refrDir);

                Vector3f reflOrig = (dotProduct(reflDir, N) < 0) ? hitPoint - scene.epsilon * N : hitPoint + scene.epsilon * N;
                Vector3f refrOrig = (dotProduct(refrDir, N) < 0) ? hitPoint - scene.epsilon * N : hitPoint + scene.epsilon * N;

                hitColor = kr * castRay(reflOrig, reflDir, scene, depth + 1)
                    + (1 - kr) * castRay(refrOrig, refrDir, scene, depth + 1);

                break;
            }
            case REFLECTION:
            {
                // TODO: 反射方向 + 递归 + fresnel 权重
                float kr = fresnel(dir, N, payload->hit_obj->ior);
                Vector3f reflDir = normalize(reflect(dir, N));
                Vector3f reflOrig = (dotProduct(reflDir, N) < 0) ? hitPoint - scene.epsilon * N : hitPoint + scene.epsilon * N;
                hitColor = kr * castRay(reflOrig, reflDir, scene, depth + 1);

                break;
            }
            default:
            {
                // TODO: 遍历所有光源，算 diffuse + specular，并用 shadow ray 判断遮挡
                Vector3f lightAmt = 0;
                Vector3f specularAmt = 0;
                Vector3f shadowOffset = (dotProduct(dir, N) < 0 ? N : -N) * scene.epsilon;
                
                for (auto& light : scene.get_lights())//遍历每一条光线
                {
                    //Lambert
                    Vector3f lightVec = light->position - hitPoint;
                    float lightDist2 = dotProduct(lightVec, lightVec);
                    Vector3f lightDir = normalize(lightVec);

                    float LdotN = std::max(0.f, dotProduct(N, lightDir));
                    if (LdotN <= 0.f) continue;

                    Vector3f shadowOrig = hitPoint + shadowOffset;

                    auto shadow_res = trace(shadowOrig, lightDir, scene.get_objects());
                    bool inShadow = shadow_res && shadow_res->tNear * shadow_res->tNear < lightDist2;

                    //Bulin-Phone
                    Vector3f viewDir = -dir;
                    Vector3f halfDir = normalize(viewDir + lightDir);
                    float p = payload->hit_obj->specularExponent;
                    float spec = powf(std::max(0.f, dotProduct(N, halfDir)), p);

                    if (!inShadow) {
                        lightAmt += light->intensity * LdotN;
                        specularAmt += light->intensity * spec;
                    }
                }
                hitColor = lightAmt * payload->hit_obj->evalDiffuseColor(st) * payload->hit_obj->Kd
                    + specularAmt * payload->hit_obj->Ks;
                break;
            }
        }
    }
    return hitColor;
}

// [comment]
// The main render function. This where we iterate over all pixels in the image, generate
// primary rays and cast these rays into the scene. The content of the framebuffer is
// saved to a file.
// [/comment]
void Renderer::Render(const Scene& scene)
{
    std::vector<Vector3f> framebuffer(scene.width * scene.height);

    float scale = std::tan(deg2rad(scene.fov * 0.5f));
    float imageAspectRatio = scene.width / (float)scene.height;

    // Use this variable as the eye position to start your rays.
    Vector3f eye_pos(0);

    int m = 0;
    for (int j = 0; j < scene.height; ++j)
    {
        for (int i = 0; i < scene.width; ++i)
        {
            // generate primary ray direction
            float x;
            float y;
            // TODO: Find the x and y positions of the current pixel to get the direction
            // vector that passes through it.
            // Also, don't forget to multiply both of them with the variable *scale*, and
            // x (horizontal) variable with the *imageAspectRatio*

			// 像素坐标以屏幕中心为原点，x轴水平向右，y轴竖直向下
			// x 从 [0, width] 映射到 [-1, 1]
            x = (2 * (i + 0.5) / (float)scene.width - 1) * imageAspectRatio * scale;
			// y 从 [0, height] 映射到 [1, -1]
            y = (1 - 2 * (j + 0.5) / (float)scene.height) * scale;

            // Don't forget to normalize this direction!
            Vector3f dir = normalize(Vector3f(x, y, -1));

            framebuffer[m++] = castRay(eye_pos, dir, scene, 0);
        }
        UpdateProgress((j + 1) / (float)scene.height);
    }

    // save framebuffer to file
    FILE* fp = fopen("binary.ppm", "wb");
    (void)fprintf(fp, "P6\n%d %d\n255\n", scene.width, scene.height);
    for (auto i = 0; i < scene.height * scene.width; ++i) {
        static unsigned char color[3];
        color[0] = (char)(255 * clamp(0, 1, framebuffer[i].x));
        color[1] = (char)(255 * clamp(0, 1, framebuffer[i].y));
        color[2] = (char)(255 * clamp(0, 1, framebuffer[i].z));
        fwrite(color, 1, 3, fp);
    }
    fclose(fp);    
}
