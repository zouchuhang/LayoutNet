function cuboidhyp_omap = generate_cuboid_from_omap(omapmore, vp, OMAP_FACTOR)

[depthorder reglabel regori rpo] = omap_depthorder(omapmore, vp, OMAP_FACTOR);
cuboidhyp_omap = omapobj2cuboidhyp(rpo, reglabel, regori, vp, OMAP_FACTOR);

