#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "draw.h"
#include <framework/trackball.h>

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableDepthOfField) {
        return;
    }

    // ...
}

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }

    // ...
}

std::array<glm::vec3, 3> constructOrthonormalBasis(const glm::vec3& r)
{
    glm::vec3 w = glm::normalize(r);
    glm::vec3 t;

    const glm::vec3 w_abs = glm::abs(w);

    if (w_abs.x < w_abs.y && w_abs.x < w_abs.z)
    {
        t = glm::vec3 { 1.0f, w.y, w.z };
    }
    else if (w_abs.y < w_abs.x && w_abs.y < w_abs.z)
    {
        t = glm::vec3 { w.x, 1.0f, w.z };
    } 
    else 
    {
        t = glm::vec3 { w.x, w.y, 1.0f };
    }
    
    glm::vec3 u = glm::normalize(glm::cross(t, w));
    glm::vec3 v = glm::cross(w, u);
    return { w, u, v };
}

std::array<float, 2> sampleDisk(RenderState& state, const float radius = 1.0f)
{
    const float r = radius * glm::sqrt(state.sampler.next_1d());
    const float angle = state.sampler.next_1d() * 2 * glm::pi<float>();
    const float x = r * glm::cos(angle);
    const float y = r * glm::sin(angle);
    return { x, y };
}


// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    const uint32_t numSamples = state.features.extra.numGlossySamples;

    if (numSamples <= 0) return;
    
    // Radius of the disk. 64 / shininess was too high for the provided scenes.
    const float radius = 0.5f / hitInfo.material.shininess;

    // "Normally" reflected ray
    Ray r = generateReflectionRay(ray, hitInfo);
    std::array<glm::vec3, 3> basis = constructOrthonormalBasis(r.direction);

    glm::vec3 accumulatedColor {};

    for (int i = 0; i < numSamples; i++)
    {
        const std::array<float, 2> point = sampleDisk(state, radius);

        Ray perturbedRay;
        perturbedRay.direction = glm::normalize(basis[0] + point[0] * basis[1] + point[1] * basis[2]);

        if (glm::dot(perturbedRay.direction, hitInfo.normal) <= 0) 
        {
            continue;
        }

        perturbedRay.origin = r.origin + 0.0001f * perturbedRay.direction;
        perturbedRay.t = std::numeric_limits<float>::max();

        accumulatedColor += hitInfo.material.ks * renderRay(state, perturbedRay, rayDepth + 1);
    }
    
    hitColor += accumulatedColor / float(numSamples);

    // Visual debug
    HitInfo hitInfoCopy = hitInfo;
    state.bvh.intersect(state, r, hitInfoCopy);
    drawRay(r, glm::vec3 { 0.5f, 0.0f, 0.8f });

    const float sphereDistFactor = 0.3f;
    const glm::vec3 sphereColor = glm::vec3 { 0.75f, 0.85f, 0.0f };
    drawSphere(r.origin + sphereDistFactor * basis[0], sphereDistFactor * radius, sphereColor);
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here
        return glm::vec3(0.f);
    } else {
        return glm::vec3(0.f);
    }
}


// TODO: Extra feature
// As an alternative to `splitPrimitivesByMedian`, use a SAH+binning splitting criterion. Refer to
// the `Data Structures` lecture for details on this metric.
// - aabb;       the axis-aligned bounding box around the given triangle set
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires splitting
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesBySAHBin(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    using Primitive = BVH::Primitive;

    return 0; // This is clearly not the solution
}