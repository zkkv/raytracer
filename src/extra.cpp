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

float perceivedLuminance(glm::vec3 colors)
{
    // Coefficients for calculation of perceived luminance
    const float c_R = 0.2126f;
    const float c_G = 0.7152f;
    const float c_B = 0.0722f;
    return c_R * colors.r + c_G * colors.g + c_B * colors.b;
}

void applyThreshold(Screen& image, const float threshold)
{
    std::vector<glm::vec3>& img = image.pixels();
    for (size_t i = 0; i < img.size(); i++)
    {
        if (perceivedLuminance(img[i]) < threshold)
        {
            img[i] = glm::vec3 { 0.0f };
        }
    }
}

void renderBloomAboveThreshold(Screen& image, const float threshold)
{
    for (int i = 0; i < image.pixels().size(); i++) 
    {
        const float luminance = perceivedLuminance(image.pixels()[i]);
        if (luminance >= threshold) 
        {
            const float mapped = (luminance - threshold) / (1 - threshold);
            // const float temp = luminance - threshold;
            // image.pixels()[i] = (1 - luminance) * glm::vec3 { 0.8f, 0.0f, 0.0f } + (luminance) * glm::vec3 { 0.85f, 0.85f, 0.85f };
            // image.pixels()[i] = (luminance) * glm::vec3 { 0.8f, 0.0f, 0.0f } + (1 - luminance) * glm::vec3 { 0.0f, 0.0f, 0.2f };
            image.pixels()[i] = (1 - mapped) * glm::vec3 { 0.8f, 0.0f, 0.0f } + (mapped)*glm::vec3 { 0.85f, 0.85f, 0.85f };
        } 
        else 
        {
            // const float temp = luminance;
            // image.pixels()[i] = (1 - temp) * glm::vec3 { 0.0f, 0.0f, 0.2f } + temp * glm::vec3 { 1.0f, 1.0f, 1.0f };
            image.pixels()[i] = glm::vec3 { 0.0f, 0.0f, 0.2f };
        }
        /*image.pixels()[i] = luminance * glm::vec3 { 0.8f, 0.0f, 0.0f } + (1.0f - luminance) * glm::vec3 {0.0f, 0.4f, 0.7f}; */
    }
}

void renderBloomBlurredMask(Screen& image, const Screen& mask, const uint32_t filterSize)
{
    std::vector<glm::vec3>& img = image.pixels();
    const std::vector<glm::vec3>& msk = mask.pixels();

    const uint32_t width = image.resolution().x;
    const uint32_t height = image.resolution().y;

    for (size_t x = 0; x < width; x++)
    {
        for (size_t y = 0; y < height; y++) 
        {
            const int index = mask.indexAt(x + filterSize, y + filterSize);
            image.setPixel(x, y, msk[index]);
        }
    }
}

Screen padBorders(const Screen& image, const uint32_t filterSize)
{
    //std::cout << "PAD START" << std::endl;
    const glm::vec3 padding { 0.0f };

    const uint32_t newWidth = image.resolution().x + 2 * filterSize;
    const uint32_t newHeight = image.resolution().y + 2 * filterSize;

    Screen padded(glm::ivec2(newWidth, newHeight), false);

    for (uint32_t x = 0; x < newWidth; x++)
    {
        for (uint32_t y = 0; y < newHeight; y++)
        {
            //std::cout << x << " " << y << " " << x - filterSize << " " << y - filterSize;
            if ((x >= filterSize) && (x < newWidth - filterSize) && (y >= filterSize) && (y < newHeight - filterSize))
            {
                //std::cout << " IN" << std::endl;
                const uint32_t index = image.indexAt(x - filterSize, y - filterSize);
                padded.setPixel(x, y, image.pixels()[index]);
            }
            else
            {
                //std::cout << " OUT" << std::endl;
                padded.setPixel(x, y, padding);
            }
        }
    }

    //std::cout << "PAD END\n" << std::endl;
    return padded;
}

uint64_t binomial(const int n, const int k)
{
    if (k == 0)
    {
        return 1;
    }
    return (n * binomial(n - 1, k - 1)) / k;
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

std::vector<glm::vec3> generateGaussianFilter(const int K)
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

    return oneDim;
}

void applyFilter1dToPixel(const size_t xOriginal, 
                        const size_t yOriginal,
                        Screen& imageCopy,
                        const std::vector<glm::vec3>& filter,
                        const size_t filterSize,
                        const bool toRow)
{
    //if (xOriginal == 799 && yOriginal == 799)
        //std::cout << "DEBUG" << std::endl;
    const size_t x = xOriginal + filterSize;
    const size_t y = yOriginal + filterSize;
    //const auto& pixelArray = image.pixels();
    std::vector<glm::vec3>& imgCopy = imageCopy.pixels();
    if (toRow)
    {

        glm::vec3 updated { 0.0f };
        for (size_t col = x - filterSize; col <= x + filterSize; col++) 
        {
            const size_t index = imageCopy.indexAt(col, y);
            updated += imgCopy[index] * filter[col - x + filterSize];
        }
        imageCopy.setPixel(x, y, updated);
    }
    else
    {
        glm::vec3 updated { 0.0f };
        for (size_t row = y - filterSize; row <= y + filterSize; row++) 
        {
            const size_t index = imageCopy.indexAt(x, row);
            updated += imgCopy[index] * filter[row - y + filterSize];
        }
        imageCopy.setPixel(x, y, updated);
    }

    //if (toRow) 
    //{

    //    glm::vec3 updated { 0.0f };
    //    for (size_t i = -filterSize; i <= filterSize; i++) {
    //        const size_t col = x + i;
    //        const size_t index = imageCopy.indexAt(col, y);
    //        updated += imgCopy[index] * filter[i + filterSize];
    //    }
    //    imageCopy.setPixel(x, y, updated);
    //} 
    //else 
    //{
    //    glm::vec3 updated { 0.0f };
    //    for (size_t i = -filterSize; i <= filterSize; i++) {
    //        const size_t row = y + i;
    //        const size_t index = imageCopy.indexAt(x, row);
    //        updated += imgCopy[index] * filter[i + filterSize];
    //    }
    //    imageCopy.setPixel(x, y, updated);
    //}

    //image.setPixel(i, j, updated);
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
    // TODO: Refactor applyFilter function
    // TODO: Remove bunch of comments
    // TODO: consider different return type for binomial()
    // TODO: Maybe replace sum with powers of 2
    // TODO: Increasing filter size makes the blur go up and to the right
    // ADVICE: Make a copy -> Set all values below threshold to zero -> 
    //  blur that image -> add two images (maybe multiply the thresheld image by a number < 1)
    // ADVICE: See render.cpp: #pragma omp parallel for schedule(guided) to parallelize the for loop (if you have time)
    // ADVICE: For visual debug show the thresheld picture
    // ADVICE: For boundaries pad with a constant color (e.g. 0)

    if (!features.extra.enableBloomEffect) 
    {
        return;
    }

    const float threshold = features.extra.bloomFilterThreshold;

    // Visual debug
    if (features.extra.enableBloomShowAboveThreshold) 
    {
        renderBloomAboveThreshold(image, threshold);
        return;
    }

    const uint32_t filterSize = features.extra.bloomFilterSize;
    //const uint32_t filterSize = 2;
    const uint32_t K = filterSize * 2 + 1;
    Screen imageCopy = padBorders(image, filterSize);
    //Screen imageCopy = image;

    // Set all values below threshold to 0
    applyThreshold(imageCopy, threshold);

    //const size_t SIZE = imageCopy.size();
    const uint32_t width = image.resolution().x;
    const uint32_t height = image.resolution().y;

    const uint32_t newWidth = width + 2 * filterSize;
    const uint32_t newHeight = height + 2 * filterSize;
    //const int variance = 1;
    //const float norm = 1 / (2 * glm::pi<float>() * variance);
    //const int K = int(2 * glm::pi<float>() * glm::sqrt(variance));
    
    //std::vector<glm::vec3> filtered;
    //filtered.reserve(SIZE);

  //  //const std::array<std::array<glm::vec3, K>, K>& filter = computeGaussianFilter<K>(K);

    std::vector<glm::vec3> filter1d = generateGaussianFilter(K);

  //  /*std::array<std::array<glm::vec3, K>, K> filter;
  //  for (int i = 0; i < K; i++)
  //  {
  //      glm::vec3 sum {0.0f};
  //      for (int j = 0; j < K; j++)
  //      {
  //          filter[i][j] = glm::vec3 { float(binomial(K - 1, j)) };
  //          sum += filter[i][j];
  //      }
  //      for (int j = 0; j < K; j++) {
  //          filter[i][j] /= sum;
  //      }
  //  }

  //  for (int i = 0; i < K; i++) 
  //  {
  //      glm::vec3 sum { 0.0f };
  //      for (int j = 0; j < K; j++) 
  //      {
  //          filter[j][i] *= glm::vec3 { float(binomial(K - 1, j)) };
  //          sum += filter[j][i];
  //      }
  //      for (int j = 0; j < K; j++) 
  //      {
  //          filter[j][i] /= sum;
  //          std::cout << filter[i][j].x << " ";
  //      }
  //      std::cout << "\n";
  //  }

  //  std::cout << std::endl;*/

  ///*  for (int i = 0; i < K; i++)
  //  {
  //      for (int j = 0; j < K; j++)
  //      {
  //          std::cout << filter[i][j].x << "\t\t"; 
  //      }
  //      std::cout << "\n";
  //  }
  //  std::cout << std::endl;*/

    for (size_t y = 0; y < height; y++)
    {
        //glm::vec3 sum {0.0f};
        for (size_t x = 0; x < width; x++)
        {
            //const int index = image.indexAt(i, j);
            //const glm::vec3& pixel = imageCopy[image.indexAt(i, j)];

            applyFilter1dToPixel(x, y, imageCopy, filter1d, filterSize, true);

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

    for (size_t x = 0; x < width; x++) 
    {
        for (size_t y = 0; y < height; y++) 
        {
            applyFilter1dToPixel(x, y, imageCopy, filter1d, filterSize, false);
        }
    }

    // Visual debug
    if (features.extra.enableBloomShowBlurredMask) {
        renderBloomBlurredMask(image, imageCopy, filterSize);
        return;
    }

    const float bloomFactor = features.extra.bloomFilterIntensity;
    for (int x = 0; x < width; x++) 
    {
        for (int y = 0; y < height; y++) 
        {
            const int origIndex = image.indexAt(x, y);
            const int padIndex = imageCopy.indexAt(x + filterSize, y + filterSize);
            const glm::vec3 stacked = image.pixels()[origIndex] + imageCopy.pixels()[padIndex] * bloomFactor;
            image.setPixel(x, y, stacked);
        }
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