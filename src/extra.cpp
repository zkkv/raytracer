#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include <iostream>

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
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...
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

float calculateAABBSurfaceArea(const AxisAlignedBox& aabb)
{
    // TODO: Test
    // TODO: make simpler
    const glm::vec3 axes = glm::abs(aabb.upper - aabb.lower);
    const float S1 = axes.x * axes.y;
    const float S2 = axes.x * axes.z;
    const float S3 = axes.y * axes.z;
    return 2 * (S1 + S2 + S3);
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
    // TODO: Test
    // TODO: Visual debugging
    // ASK: What should happen if nPrimitives < nBins?

    /* DEBUG START */
    /*const size_t DEBUG_SIZE = 17;
    int arr[DEBUG_SIZE];
    for (int i = 0; i < DEBUG_SIZE; i++)
    {
        arr[i] = i;
    }
    std::span<int, DEBUG_SIZE> p = std::span(arr);
    const size_t N = DEBUG_SIZE;*/
    /* DEBUG END */

    using Primitive = BVH::Primitive;

    const size_t N = primitives.size();
    size_t nBins = 10;

    // Not sure what to do in this case
    if (N < nBins)
    {
        nBins = N;
        //return splitPrimitivesByMedian(aabb, axis, primitives);
    }

    // Sort
    splitPrimitivesByMedian(aabb, axis, primitives);

    // The last bin can be bigger because primitives.size() might be not divisible by nBins.
    const size_t binSize = N / nBins;
    const size_t lastBinSize = binSize + N % nBins;

    // We actually don't care about this since this is constant.
    //const float outerSurfaceArea = calculateAABBSurfaceArea(aabb);

    float minCost = std::numeric_limits<float>::max();
    size_t minIndex = 1;

    for (size_t i = 1; i < nBins; i++)
    {
        const size_t iSplit = i * binSize;
        const size_t leftSize = iSplit;
        const size_t rightSize = N - leftSize;
        const float leftSurfaceArea = calculateAABBSurfaceArea(computeSpanAABB(primitives.subspan(0, leftSize)));
        const float rightSurfaceArea = calculateAABBSurfaceArea(computeSpanAABB(primitives.subspan(iSplit, rightSize)));
        const float cost = leftSurfaceArea * leftSize + rightSurfaceArea * rightSize;
        if (cost < minCost)
        {
            minCost = cost;
            minIndex = i;
        }

        /* DEBUG START */
        /*{
            std::cout << "i = " << i << ", iSplit = " << iSplit << ", leftSize = " << leftSize << ", rightSize = " << rightSize;
            std::cout << "\nLEFT:\n";

            const auto leftSpan = p.subspan(0, leftSize);
            const auto rightSpan = p.subspan(iSplit, rightSize);

            for (const auto& elem : leftSpan) {
                std::cout << elem << " ";
            }
            std::cout << "\nRIGHT:\n";
            for (const auto& elem : rightSpan) {
                std::cout << elem << " ";
            }
            std::cout << "\n\n";
        }*/
        /* DEBUG END */
    }

    return minIndex * binSize;

    //return 0; // This is clearly not the solution
}