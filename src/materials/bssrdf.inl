/*
	3/15/2026

	This is only a trivial override.
	The real implementation is in bssrdf.cpp

	I do it this way because the path_tracing.cpp pipeline doens't fit bssrdf sampling and eval etc.
	So I create a new branch in path_tracing.cpp to deal with bssrdf mtl.

*/

Spectrum eval_op::operator()(const BSSRDF& bssrdf) const {
	//Spectrum L_single_scatter = L_single();
	//Spectrum L_multi_scatter = L_diffusion(scene, vertex, bssrdf, dir_in, init_pcg32(0,0));

	//return L_single_scatter + L_multi_scatter;
	return make_const_spectrum(0);
}

Real pdf_sample_bsdf_op::operator()(const BSSRDF& bssrdf) const {
	return 0;
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const BSSRDF& bssrdf) const {
	return BSDFSampleRecord{ Vector3(0,0,0), 1, 1 };
}

TextureSpectrum get_texture_op::operator()(const BSSRDF& bssrdf) const {
	return bssrdf.base_color;
}



