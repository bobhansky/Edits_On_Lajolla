#include "material.h"
#include "intersection.h"

inline Vector3 sample_cos_hemisphere(const Vector2 &rnd_param) {
    Real phi = c_TWOPI * rnd_param[0];
    Real tmp = sqrt(std::clamp(1 - rnd_param[1], Real(0), Real(1)));
    return Vector3{
        cos(phi) * tmp, sin(phi) * tmp,
        sqrt(std::clamp(rnd_param[1], Real(0), Real(1)))
    };
}



#include "materials/lambertian.inl"
#include "materials/roughplastic.inl"
#include "materials/roughdielectric.inl"
#include "materials/disney_diffuse.inl"
#include "materials/disney_metal.inl"
#include "materials/disney_glass.inl"
#include "materials/disney_clearcoat.inl"
#include "materials/disney_sheen.inl"
#include "materials/disney_bsdf.inl"
#include "materials/bssrdf.inl"

Spectrum eval(const Material &material,
              const Vector3 &dir_in,
              const Vector3 &dir_out,
              const PathVertex &vertex,
              const TexturePool &texture_pool,
              TransportDirection dir) {
    return std::visit(eval_op{dir_in, dir_out, vertex, texture_pool, dir}, material);
}

std::optional<BSDFSampleRecord>
sample_bsdf(const Material &material,
            const Vector3 &dir_in,
            const PathVertex &vertex,
            const TexturePool &texture_pool,
            const Vector2 &rnd_param_uv,
            const Real &rnd_param_w,
            TransportDirection dir) {
    return std::visit(sample_bsdf_op{
        dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w, dir}, material);
}

Real pdf_sample_bsdf(const Material &material,
                     const Vector3 &dir_in,
                     const Vector3 &dir_out,
                     const PathVertex &vertex,
                     const TexturePool &texture_pool,
                     TransportDirection dir) {
    return std::visit(pdf_sample_bsdf_op{
        dir_in, dir_out, vertex, texture_pool, dir}, material);
}

TextureSpectrum get_texture(const Material &material) {
    return std::visit(get_texture_op{}, material);
}


Texture<Spectrum> get_normalMap(const Material& material) {
    return std::visit(get_normalMap_op{}, material);
}

bool has_normal(const Material& material) {
    return std::visit(has_normal_op{}, material);
}