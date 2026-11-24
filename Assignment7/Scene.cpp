//
// Created by G枚ksu G眉vendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // ===== 1. 这条光线打到了什么 =====
    //   未命中     → 返回黑色
    //   命中光源   → 只有“相机直接看到光源”才返回它的 emission；
    //                弹射之后打到光源的情况已经由 L_dir 计入，这里不能再算一次
    //   命中表面   → 往下算 L_dir + L_indir
    Intersection inter = intersect(ray);
    if (!inter.happened) {
        return Vector3f(0.f);
    }
    else if (inter.m->hasEmission()) {
        // TODO: 这里该返回什么？用什么条件区分“相机直接看到”和“弹射后打到”？
        if (depth == 0) {
            return Vector3f(inter.m->getEmission());
        }
        else {
            return Vector3f(0.f);
        }
    }

    // ===== 2. 直接光 L_dir =====
    //   输入量：p = inter.coords，N = inter.normal
    //           xs / pdf_light ← 由 sampleLight(light_pos, pdf_light) 给出
    //   遮挡判断：从 p 朝 xs 打一条阴影射线，最近命中距离 ≈ |xs - p| 才算没被挡
    //   权重：f_r * cos(theta_p) * cos(theta_s) / (|xs - p|^2 * pdf_light)
    //   坑：Material::eval / pdf 的第二个参数必须是入射方向，传反会静默返回 0
    Vector3f p = inter.coords;
    Vector3f n = inter.normal;

    float pdf_light;
    Intersection light_pos;
    sampleLight(light_pos, pdf_light);

    Vector3f ws = normalize(light_pos.coords - p);
    float dist = (p - light_pos.coords).norm();

    float eps= 1e-3f;
    Vector3f OffsetDir = (dotProduct(n, ws) > 0) ? n : -n;
    Vector3f ShadowOrigin = p + eps * OffsetDir;
    Ray shadowRay(ShadowOrigin, ws);
    
    Intersection s = intersect(shadowRay);
    Vector3f L_dir(0.f);

    // 遮挡判断：s.distance 是阴影射线的最近命中距离
    // 没遮挡时 s.distance ≈ dist（差一个 eps），Intersection 默认 distance = DBL_MAX
    bool blocked = (dist - s.distance) > eps;

    if (!blocked) {
        float cos_p = std::max(0.f, dotProduct(n, ws));                   // 表面法线 · 指向光源
        float cos_s = std::max
            (dotProduct(light_pos.normal, -ws),0.f);    // 光源法线 · 指向表面

        L_dir = light_pos.emit                               // L_e(x_s)
            * inter.m->eval(ws, -ray.direction, n)         // f_r
            * cos_p * cos_s                                // 两个余弦
            / (dist * dist * pdf_light);                   // 面积→立体角的雅可比 × pdf
    }
    
    // ===== 3. 间接光 L_indir =====
    //   方向：Material::sample 采样得到，不是自己解算
    //   权重：f_r * cos(theta_i) / pdf(wi) / P_RR
    //   轮盘：以 P_RR 的概率继续，结果除以 P_RR；它不是深度计数器
    //   递归命中光源时跳过，避免和 L_dir 重复计入
    Vector3f L_indir(0.f);

    if (get_random_float() < RussianRoulette) {
        // 1. 采样方向：用材质，不自己解算
        Vector3f wi = inter.m->sample(-ray.direction, n);

        // 2. 这个方向的 pdf
        float pdf_i = inter.m->pdf(wi, -ray.direction, n);

        if (pdf_i > 1e-6f) {
            // 3. 沿 wi 的半球偏移，避免自相交
            Vector3f offsetDir2 = (dotProduct(n, wi) > 0) ? n : -n;

            // 4. 递归
            Vector3f Li = castRay(Ray(p + eps * offsetDir2, wi), depth + 1);

            // 5. MC 权重
            L_indir = inter.m->eval(wi, -ray.direction, n)   // f_r
                * dotProduct(n, wi)                       // cos(theta_i)
                * Li                                      // L_i
                / pdf_i                                   // 1 / pdf
                / RussianRoulette;                  // 1 / P_RR
        }
    }

    // ===== 4. 返回 =====
    return L_dir + L_indir;
}