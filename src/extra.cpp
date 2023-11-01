#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "draw.h"
#include <framework/trackball.h>
#include "texture.h"
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
    // Renders the pixels which are above the threshold in red-white gradient.
    // All other values are set to dark blue.
    std::vector<glm::vec3>& img = image.pixels();
    for (size_t i = 0; i < img.size(); i++) 
    {
        const float luminance = perceivedLuminance(img[i]);
        if (luminance >= threshold) 
        {
            const float mapped = (luminance - threshold) / (1 - threshold);
            img[i] = (1 - mapped) * glm::vec3 { 0.8f, 0.0f, 0.0f } + (mapped) * glm::vec3 { 0.85f, 0.85f, 0.85f };
        } 
        else 
        {
            img[i] = glm::vec3 { 0.0f, 0.0f, 0.2f };
        }
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
            const int indexPadded = mask.indexAt(x + filterSize, y + filterSize);
            image.setPixel(x, y, msk[indexPadded]);
        }
    }
}

Screen padBorders(const Screen& image, const uint32_t filterSize)
{
    // Padding borders with black color (default color in Screen constructor)
    const uint32_t newWidth = image.resolution().x + 2 * filterSize;
    const uint32_t newHeight = image.resolution().y + 2 * filterSize;

    Screen padded(glm::ivec2(newWidth, newHeight), false);

    for (size_t x = 0; x < newWidth; x++)
    {
        for (size_t y = 0; y < newHeight; y++)
        {
            if ((x >= filterSize) && (x < newWidth - filterSize) && (y >= filterSize) && (y < newHeight - filterSize))
            {
                const size_t index = image.indexAt(x - filterSize, y - filterSize);
                padded.setPixel(x, y, image.pixels()[index]);
            }
        }
    }
    return padded;
}

uint64_t binomial(const uint32_t n, const uint32_t k)
{
    uint32_t k_min = k > n - k ? n - k : k;
    uint64_t coef = 1;

    // Multiplicative formula
    for (uint32_t i = 1; i <= k_min; i++) 
    {
        coef = (coef * (n - i + 1)) / i;
    }

    return coef; 
}

std::vector<glm::vec3> generateGaussianFilter(const uint32_t filterSize)
{
    const uint32_t radius = filterSize * 2 + 1;
    std::vector<glm::vec3> filter;
    filter.reserve(radius);
    glm::vec3 sum { 0.0f };

    for (size_t i = 0; i < radius; i++)
    {
        const glm::vec3 computed = glm::vec3 { float(binomial(radius - 1, i)) };
        filter.push_back(computed);
        sum += computed;
    }

    for (size_t i = 0; i < radius; i++) 
    {
        filter[i] /= sum;
    }

    return filter;
}

void applyFilter1dToPixel(const size_t xOriginal, 
                        const size_t yOriginal,
                        const Screen& imagePadded,
                        Screen& imageCopy,
                        const std::vector<glm::vec3>& filter,
                        const uint32_t filterSize,
                        const bool toRow)
{
    const size_t x = xOriginal + filterSize;
    const size_t y = yOriginal + filterSize;
    const std::vector<glm::vec3>& imgPadded = imagePadded.pixels();

    glm::vec3 updated { 0.0f };
    if (toRow)
    {
        for (size_t col = x - filterSize; col <= x + filterSize; col++) 
        {
            const size_t index = imagePadded.indexAt(col, y);
            updated += imgPadded[index] * filter[col - x + filterSize];
        }
    }
    else
    {
        for (size_t row = y - filterSize; row <= y + filterSize; row++) 
        {
            const size_t index = imagePadded.indexAt(x, row);
            updated += imgPadded[index] * filter[row - y + filterSize];
        }
    }
    imageCopy.setPixel(x, y, updated);
}

// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    /*
        SOURCES:
        Perceived luminance: https://stackoverflow.com/a/596243/15236567
        
        Multiplicative formula for binomial coefficients: https://stackoverflow.com/a/15302394/15236567
        
        Convolution filters: Marschner, S.; Shirley, P. Fundamentals of Computer Graphics, Fourth.; 
        CRC Press, Taylor & Francis Group: Boca Raton, FL, 2015, p. 190
        
        Binomial filters: http://www.cse.yorku.ca/~kosta/CompVis_Notes/binomial_filters.pdf.old
    */

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

    // This value MUST be limited by 31, otherwise integer overflow happens. I limit it to 30.
    const uint32_t filterSize = features.extra.bloomFilterSize;

    // imagePadded is used as an untouched copy in filtering
    Screen imagePadded = padBorders(image, filterSize);

    // Set all values below threshold to 0
    applyThreshold(imagePadded, threshold);

    // imageCopy is the one that changes
    Screen imageCopy = imagePadded;
    
    // One-dimensional binomial (Gaussian) filter
    std::vector<glm::vec3> filter1d = generateGaussianFilter(filterSize);

    const uint32_t width = image.resolution().x;
    const uint32_t height = image.resolution().y;

    // Filter is separable, so we can first go through all rows, then all columns
    for (size_t y = 0; y < height; y++)
    {
        for (size_t x = 0; x < width; x++)
        {
            applyFilter1dToPixel(x, y, imagePadded, imageCopy, filter1d, filterSize, true);
        }
    }

    // Second traversal should be done on the result of the first one
    imagePadded = imageCopy;

    for (size_t x = 0; x < width; x++) 
    {
        for (size_t y = 0; y < height; y++) 
        {
            applyFilter1dToPixel(x, y, imagePadded, imageCopy, filter1d, filterSize, false);
        }
    }

    // Visual debug
    if (features.extra.enableBloomShowBlurredMask) 
    {
        renderBloomBlurredMask(image, imageCopy, filterSize);
        return;
    }

    const float bloomFactor = features.extra.bloomFilterIntensity;
    for (int x = 0; x < width; x++) 
    {
        for (int y = 0; y < height; y++) 
        {
            const int indexOriginal = image.indexAt(x, y);
            const int indexPadded = imageCopy.indexAt(x + filterSize, y + filterSize);
            const glm::vec3 combined = image.pixels()[indexOriginal] + imageCopy.pixels()[indexPadded] * bloomFactor;
            image.setPixel(x, y, combined);
        }
    }

    //image.writeBitmapToFile("D:/final_image.bmp");
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
    /*
        SOURCES:
        Sampling uniformly on a disk: https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly
       
        Creating perturbed rays: Marschner, S.; Shirley, P. Fundamentals of Computer Graphics, Fourth.; 
        CRC Press, Taylor & Francis Group: Boca Raton, FL, 2015, chapter 13.4.4.
        
        Constructing an orthonormal basis from a single vector: Marschner, S.; Shirley, P. Fundamentals of Computer Graphics, Fourth.; 
        CRC Press, Taylor & Francis Group: Boca Raton, FL, 2015, chapter 2.4.6.
    */
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
        float x = std::fabs(ray.direction.x), y = std::fabs(ray.direction.y), z = std::fabs(ray.direction.z);
        float maxComponent = std::max(x, std::max(y, z));
        glm::vec3 r = ray.direction / maxComponent; // [-1, 1]
        glm::vec3 coords = (r + glm::vec3(1.0f)) / 2.0f; // [0, 1]

        // +- 1 => choose face
        // take other 2 coords and sample from face
        float one = 1.0f - FLT_EPSILON;
        float u, v;
        
        /* texture is 4 squares wide and 3 squares:

                   UP
            LEFT FRONT RIGHT BACK
                  DOWN
        */
        // Some 'coords' components need to be flipped to account for the way a cube is unfolded onto a flat plane (some faces are inverted)
        if (r.x > one) { // right
            u = coords.z / 4.0f + 2.0f / 4.0f;
            v = coords.y / 3.0f + 1.0f / 3.0f;
        } else if (r.x < -one) { // left
            u = (1 - coords.z) / 4.0f;
            v = coords.y / 3.0f + 1.0f / 3.0f;
        } else if (r.y > one) { // up
            u = coords.x / 4.0f + 1.0f / 4.0f;
            v = coords.z / 3.0f + 2.0f / 3.0f;
        } else if (r.y < -one) { // down
            u = coords.x / 4.0f + 1.0f / 4.0f;
            v = (1 - coords.z) / 3.0f;
        } else if (r.z < -one) { // front
            u = coords.x / 4.0f + 1.0f / 4.0f;
            v = coords.y / 3.0f + 1.0f / 3.0f;
        } else /*if (r.z > one)*/ { // back
            u = (1 - coords.x) / 4.0f + 3.0f / 4.0f;
            v = coords.y / 3.0f + 1.0f / 3.0f;
        }
        if (u < 0.0f || u > 1.0f || v < 0.0f || v > 1.0f)
            std::cout << "OUT OF BOUNDS\n";
        glm::vec2 mapTexCoords = glm::vec2(u, v);

        if (state.features.enableBilinearTextureFiltering) 
            return sampleTextureBilinear(state.scene.environmentMap, mapTexCoords);

        return sampleTextureNearest(state.scene.environmentMap, mapTexCoords);

    } else {
        return glm::vec3(0.f);
    }
}

size_t determineBucketIndex(const size_t i, const size_t nBins, uint32_t axis, std::span<BVH::Primitive> primitives)
{
    const glm::vec3 subintervalLen = computePrimitiveCentroid(primitives[i]) - computePrimitiveCentroid(primitives[0]);
    const glm::vec3 intervalLen = computePrimitiveCentroid(primitives[primitives.size() - 1]) - computePrimitiveCentroid(primitives[0]);

    size_t idx;
    if (axis == 0)
    {
        idx = nBins * (subintervalLen.x) / (intervalLen.x);
    }
    else if (axis == 1)
    {
        idx = nBins*(subintervalLen.y) / (intervalLen.y);
    } 
    else 
    {
        idx = nBins*(subintervalLen.z) / (intervalLen.z);
    }

    // Clamp the index to not get out of bounds
    if (idx == nBins) 
    {
        idx--;
    }
    return idx;
}

float calculateAABBSurfaceArea(const AxisAlignedBox& aabb)
{
    const glm::vec3 axes = glm::abs(aabb.upper - aabb.lower);
    return 2 * (axes.x * axes.y + axes.x * axes.z + axes.y * axes.z);
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
    /*
       SOURCES:
       Surface area heuristic with binning: M. Pharr, J. Wenzel, and G. Humphreys. Physically Based Rendering, Second Edition: 
       From Theory To Implementation. Morgan Kaufmann Publishers Inc., 2nd edition, chapter 4.4.2.

       General information: TU Delft Computer Graphics course, lecture 9.
    */
    using Primitive = BVH::Primitive;

    const size_t N = primitives.size();
    size_t nBins = 50;

    if (N < nBins)
    {
        nBins = N;
    }

    // Sort
    splitPrimitivesByMedian(aabb, axis, primitives);

    std::vector<BVH::Bin> bins;
    bins.resize(nBins);

    // Assign each primitive to a bin
    for (size_t i = 0; i < N; i++)
    {
        const size_t bIndex = determineBucketIndex(i, nBins, axis, primitives);
        bins[bIndex].binPrimitives.push_back(primitives[i]);
    }

    // Determine where each bin starts
    size_t relStart = 0;
    for (size_t i = 0; i < nBins; i++)
    {
        bins[i].start = relStart;
        relStart += bins[i].binPrimitives.size();
    }

    const float outerSurfaceArea = calculateAABBSurfaceArea(aabb);

    float minCost = std::numeric_limits<float>::max();
    size_t minIndex = -1;
    
    // Intersection cost is assumed to be 1 for all calculations below
    const float baseCost = N;

    // Traversal cost can be tuned to get the best results
    const float travCost = 1.5f;
    
    // Finding min cost by considering splits after each bucket
    for (size_t i = 0; i < nBins - 1; i++) 
    {
        // Left subdivision
        const size_t leftSize = bins[i + 1].start;
        const float leftSurfaceArea = calculateAABBSurfaceArea(computeSpanAABB(primitives.subspan(0, leftSize)));

        // Right subdivision
        const size_t rightSize = N - leftSize;
        const float rightSurfaceArea = calculateAABBSurfaceArea(computeSpanAABB(primitives.subspan(leftSize, rightSize)));

        const float cost = travCost + (leftSurfaceArea * leftSize + rightSurfaceArea * rightSize) / outerSurfaceArea;
        if (cost < minCost) 
        {
            minCost = cost;
            minIndex = leftSize;
        }
    }

    if (minCost >= baseCost)
    {
        return -1;
    }

    return minIndex;
}