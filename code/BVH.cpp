#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    root = recursiveBuild(primitives);

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    printf(
        "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
        hrs, mins, secs);
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        return node;
    }
    else {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim) {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x <
                       f2->getBounds().Centroid().x;
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y <
                       f2->getBounds().Centroid().y;
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z <
                       f2->getBounds().Centroid().z;
            });
            break;
        }

        auto beginning = objects.begin();
        auto middling = objects.begin() + (objects.size() / 2);
        auto ending = objects.end();

        auto leftshapes = std::vector<Object*>(beginning, middling);
        auto rightshapes = std::vector<Object*>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
    }

    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    std::array<int, 3> dirIsNeg = {ray.direction.x > 0, ray.direction.y > 0, ray.direction.z > 0};

    if (!node->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg)) {
        return Intersection();
    }

    if (node->left == nullptr && node->right == nullptr) {
        return node->object->getIntersection(ray);
    }

    // Intersection hit1 = getIntersection(node->left, ray);
    // Intersection hit2 = getIntersection(node->right, ray);

    // return hit1.distance < hit2.distance ? hit1 : hit2;

    // 计算光线进入左右子节点包围盒的距离，用于决定遍历顺序和剪枝
    float t_enter_left  = std::numeric_limits<float>::max();
    float t_enter_right = std::numeric_limits<float>::max();

    if (node->left && node->left->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg)) {
        // 快速估算 entry distance：取光线原点到包围盒中心的距离
        t_enter_left = dotProduct(node->left->bounds.Centroid() - ray.origin, ray.direction);
    }
    if (node->right && node->right->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg)) {
        t_enter_right = dotProduct(node->right->bounds.Centroid() - ray.origin, ray.direction);
    }

    // 先遍历距离近的子节点（front-to-back 遍历）
    BVHBuildNode* first  = (t_enter_left <= t_enter_right) ? node->left  : node->right;
    BVHBuildNode* second = (t_enter_left <= t_enter_right) ? node->right : node->left;
    float t_second = std::min(t_enter_left, t_enter_right);  // 较远子节点的进入距离

    Intersection hit = getIntersection(first, ray);

    // 剪枝：near 子节点找到的交点距离 < far 子节点的进入距离 → 跳过 far
    if (hit.happened && hit.distance < t_second) {
        return hit;
    }

    Intersection hit2 = getIntersection(second, ray);
    if (hit2.happened && hit2.distance < hit.distance) {
        return hit2;
    }

    return hit;
}