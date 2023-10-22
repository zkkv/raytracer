#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include <framework/trackball.h>
#include <iostream>

void printVector(const glm::vec3& vector, const std::string& str = "")
{
    std::cout << str << "(" << vector.x << ", " << vector.y << ", " << vector.z << ")" << std::endl;
}

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

uint64_t binomial(const int n, const int k)
{
    if (k == 0)
    {
        return 1;
    }
    return (n * binomial(n - 1, k - 1)) / k;
}

float perceivedLuminance(glm::vec3 colors)
{
    // Coefficients for calculation of perceived luminance
    const float c_R = 0.2126f;
    const float c_G = 0.7152f;
    const float c_B = 0.0722f;
    return c_R * colors.r + c_G * colors.g + c_B * colors.b;
}
//
//template<int K>
//std::array<std::array<glm::vec3, K>, K>& computeGaussianFilter(int K)
//{
//    std::array<std::array<glm::vec3, K>, K> filter;
//    for (int i = 0; i < K; i++)
//    {
//        for (int j = 0; k < K; j++)
//        {
//
//        }
//    }
//    return filter;
//}

std::vector<std::vector<glm::vec3>> outerProduct(std::vector<glm::vec3>& vector)
{
    std::vector<std::vector<glm::vec3>> result;
    const size_t SIZE = vector.size();
    result.resize(SIZE);

    for (int i = 0; i < SIZE; i++)
    {
        result[i].resize(SIZE);
    }

    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
        {
            result[i][j] = vector[i] * vector[j];
        }
    }
    return result;
}

std::vector<std::vector<glm::vec3>> generateGaussianFilter(const int K)
{
    std::vector<glm::vec3> oneDim;
    oneDim.reserve(K);
    glm::vec3 sum { 0.0f };
    for (int i = 0; i < K; i++)
    {
        const glm::vec3 computed = glm::vec3 { float(binomial(K - 1, i)) };
        oneDim.push_back(computed);
        sum += computed;
    }

    for (int i = 0; i < K; i++) 
    {
        oneDim[i] /= sum;
    }

    return outerProduct(oneDim);
}

void applyFilterToPixel(const size_t i, 
                        const size_t j, 
                        Screen& image, 
                        const std::vector<glm::vec3>& imageCopy,
                        const std::vector<std::vector<glm::vec3>>& filter,
                        const size_t filterSize)
{
    //const auto& pixelArray = image.pixels();
    glm::vec3 updated { 0.0f };
    for (int y = i - filterSize; y <= i + filterSize; y++)
    {
        for (int x = j - filterSize; x <= j + filterSize; x++)
        {
            const size_t index = image.indexAt(y, x);
            updated += imageCopy[index] * filter[i - y + filterSize][j - x + filterSize];
        }
    }
    image.setPixel(i, j, updated);
}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    // TODO: Test
    // TODO: visual debug
    // TODO: thresholds
    // TODO: Handle boundaries
    // TODO: Make separable

    if (!features.extra.enableBloomEffect) {
        return;
    }

    const std::vector<glm::vec3> imageCopy = image.pixels();
    const size_t SIZE = imageCopy.size();
    const int width = image.resolution().x;
    const int height = image.resolution().y;
    //const int variance = 1;
    //const float norm = 1 / (2 * glm::pi<float>() * variance);
    //const int K = int(2 * glm::pi<float>() * glm::sqrt(variance));
    const float threshold = 0.7f;
    const int filterSize = 5;
    const int K = filterSize * 2 + 1;

    std::vector<glm::vec3> filtered;
    filtered.reserve(SIZE);

    //const std::array<std::array<glm::vec3, K>, K>& filter = computeGaussianFilter<K>(K);

    std::vector<std::vector<glm::vec3>> filter = generateGaussianFilter(K);

    /*std::array<std::array<glm::vec3, K>, K> filter;
    for (int i = 0; i < K; i++)
    {
        glm::vec3 sum {0.0f};
        for (int j = 0; j < K; j++)
        {
            filter[i][j] = glm::vec3 { float(binomial(K - 1, j)) };
            sum += filter[i][j];
        }
        for (int j = 0; j < K; j++) {
            filter[i][j] /= sum;
        }
    }

    for (int i = 0; i < K; i++) 
    {
        glm::vec3 sum { 0.0f };
        for (int j = 0; j < K; j++) 
        {
            filter[j][i] *= glm::vec3 { float(binomial(K - 1, j)) };
            sum += filter[j][i];
        }
        for (int j = 0; j < K; j++) 
        {
            filter[j][i] /= sum;
            std::cout << filter[i][j].x << " ";
        }
        std::cout << "\n";
    }

    std::cout << std::endl;*/

  /*  for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < K; j++)
        {
            std::cout << filter[i][j].x << "\t\t"; 
        }
        std::cout << "\n";
    }
    std::cout << std::endl;*/

    for (int i = filterSize; i < height - filterSize; i++)
    {
        //glm::vec3 sum {0.0f};
        for (int j = filterSize; j < width - filterSize; j++)
        {
            const int index = image.indexAt(i, j);
            const glm::vec3& pixel = imageCopy[index];
            if (perceivedLuminance(pixel) < threshold)
            {
                continue;
            }

            applyFilterToPixel(i, j, image, imageCopy, filter, filterSize);

            //printVector(imageCopy[index]);
            //printVector(image.pixels()[index]);
            //std::cout << std::endl;
            
            //const float factor = binomial(width, i);

            //const float debug = glm::length(pixel);

            //std::cout << debug << std::endl;

            //const glm::vec3 newVal = pixel * factor;

            //image.setPixel(i, j, newVal);
            //sum += newVal;
        }
        
         //for (int j = 0; j < width; j++) 
         //{
         //   const int index = image.indexAt(i, j);
         //   image.setPixel(i, j, pixelsArray[index] / sum);
         //}
    }
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