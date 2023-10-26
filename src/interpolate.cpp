#include "interpolate.h"
#include <glm/geometric.hpp>

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    glm::vec3 first = v1 - v0;
    glm::vec3 second = v2 - v0;
    glm::vec3 point = p - v0;

    float beta = (glm::dot(second, second) * glm::dot(first, point) - glm::dot(first, second) * glm::dot(second, point)) / (glm::dot(first, first) * glm::dot(second, second) - glm::dot(first, second) * glm::dot(first, second));
    float gamma = (glm::dot(first, first) * glm::dot(second, point) - glm::dot(first, second) * glm::dot(first, point)) / (glm::dot(first, first) * glm::dot(second, second) - glm::dot(first, second) * glm::dot(first, second));

    return glm::vec3 { 1.f - beta - gamma, beta, gamma };
}

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    // TODO: implement this function.
    return bc[0] * n0 + bc[1] * n1 + bc[2] * n2;
}

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
// TODO: implement this function.
    return bc[0] * t0 + bc[1] * t1 + bc[2] * t2;
}
