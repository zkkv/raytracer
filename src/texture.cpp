#include "texture.h"
#include "render.h"
#include <framework/image.h>

glm::vec3 getPixelFromCoordinates(const Image& image, const int xCoord, const int yCoord) 
{
    // Extracts the pixel at (xCoord, yCoord) from the 1D array of the image

    return image.pixels[image.pixels.size() - (yCoord * image.width + xCoord) - 1];
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
    int xCoord = glm::round(texCoord[0] * (image.width - 1));
    int yCoord = glm::round(texCoord[1] * (image.height - 1));
    
    return getPixelFromCoordinates(image, xCoord, yCoord);
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)


    // Points to interpolate from: TL = top left, TR = top right, BL = bottom left, BR = bottom right
    int xCoordTL = glm::floor(texCoord[0] * (image.width - 1));
    int yCoordTL = glm::ceil(texCoord[1] * (image.height - 1));

    int xCoordTR = glm::ceil(texCoord[0] * (image.width - 1));
    int yCoordTR = glm::ceil(texCoord[1] * (image.height - 1));

    int xCoordBL = glm::floor(texCoord[0] * (image.width - 1));
    int yCoordBL = glm::floor(texCoord[1] * (image.height - 1));

    int xCoordBR = glm::ceil(texCoord[0] * (image.width - 1));
    int yCoordBR = glm::floor(texCoord[1] * (image.height - 1));

    // Coordinates of the point we interpolate for
    float xCoord = texCoord[0] * (image.width - 1);
    float yCoord = texCoord[1] * (image.height - 1);
    
    float alpha = xCoord - xCoordBL;
    float beta = yCoord - yCoordBL;

    // Compute areas
    float coefTL = (1.f - alpha) * beta;
    float coefTR = alpha * beta;
    float coefBL = (1.f - alpha) * (1.f - beta);
    float coefBR = alpha * (1.f - beta);

    // Compute color contributed by each corner pixel
    glm::vec3 result = coefTL * getPixelFromCoordinates(image, xCoordTL, yCoordTL);
    result += coefTR * getPixelFromCoordinates(image, xCoordTR, yCoordTR);
    result += coefBL * getPixelFromCoordinates(image, xCoordBL, yCoordBL);
    result += coefBR * getPixelFromCoordinates(image, xCoordBR, yCoordBR);

    return result;
}