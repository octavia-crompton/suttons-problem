Baldocchi et al. (2016) show annual evaporation was highest when the flooded rice field was small and isolated, and decreased as flooded area expanded (interpreted as diminishing oasis/advective enhancement at larger spatial scale).

RAPID experiment (Idaho alfalfa): regional advection with negative sensible heat flux — de Bruin et al. (2005) report oasis signatures over extensive irrigated alfalfa in an arid environment (daytime negative H and LE exceeding net radiation on some days).

Crau field experiment (France): dry upwind → irrigated grass downwind, strong local advection — Kroon & de Bruin (1995) is a foundational field dataset on turbulence/scalar behavior under strong local advection across a sharp surface transition.

North China Plains irrigated wheat: regional + local advection and eddy-diffusivity inequality — Lee et al. (2004) revisits how advection alters the relationship between heat and water-vapor transport (and what that implies for ET estimation by different methods).

Advective inversion over paddy: turbulence + spectra under “disturbed” surface layer — McNaughton & Laubach (2000) analyze wind/scalar spectra near the base of an advective inversion formed over a paddy downwind of a dry region (helpful for mechanistic arguments about when “advection matters” vs. when it’s mainly an outer-layer disturbance).

Energy-balance closure + advection across sites (ADVEX) — Moderow et al. (2021) is a useful multi-site reference: advective terms can be measurable, but impacts on closure are site-specific (good for the “sometimes yes, sometimes no” section).

Orchard oasis-advection example (irrigated desert orchard) — Stoughton et al. (2002) (pecan orchard, Las Cruces NM) explicitly discusses regional (“oasis”) advection signatures and turbulence/flux-profile behavior.

Large eddies and flux convergence/divergence over irrigated cotton (EBEX context) — **Gao et al. (2016)** is a strong companion to EBEX-style discussions: large eddies can modulate apparent convergence/divergence in latent heat flux in a “disturbed” surface layer.



- The cotton field had very high evapotranspiration because of its moist soil following patch-by-patch irrigation, and this created spatial heterogeneity in ET across the field. That soil moisture heterogeneity was actually part of what generated the large eddies in the first place — a finding inherited from earlier work on the same dataset (Zhang et al., 2010).
- The latent heat flux (LE) they analyze is essentially the energy expression of ET, so in that sense the entire paper is about ET transport — just framed in energy terms rather than water volume terms.
Overall, both advection and ET heterogeneity are treated more as background context that sets the scene, while large eddies are positioned as the primary mechanism driving the flux differences between heights.


Surface resistance + advection (theory) — **Philip (1987)** directly tackles how advection interacts with evaporation and surface resistance (nice for linking your “leading-edge enhancement” framing to resistance feedbacks).

Modeling/measurement approaches to “missing transport” — **Mauder et al. (2021)** (LES-based correction options for large-scale transport not captured by single towers) is good background when you discuss why leading-edge signals are easy to miss (or misattribute)

- When scientists use eddy-covariance towers to measure how much heat and water vapor move between the land surface and the atmosphere, they consistently get answers that are too low — typically by 10–30%. The turbulent heat fluxes measured at the tower are smaller than the energy available at the surface (net radiation minus ground heat flux). This is called the "energy balance closure problem" and has been an unresolved headache in micrometeorology for decades.

- The main culprit identified in this paper is something called dispersive fluxes, which arise from large, organized circulation patterns in the atmosphere — think of slowly rotating convective cells and rolls that span the entire depth of the atmospheric boundary layer (roughly 1–2 km on a sunny day). These large structures transport heat and moisture vertically just like smaller turbulence does, but because a single tower only samples one point in space, it completely misses them. The tower captures the local turbulent eddies just fine but is blind to these kilometer-scale organized motions.

- An earlier study (De Roo et al., 2018) used large-eddy simulations — very detailed computer models of the lower atmosphere — to figure out how large these missed dispersive fluxes are under different atmospheric conditions. The result was a set of correction equations that take in two key variables: how strong the surface friction is relative to the buoyancy-driven convection (expressed as u*/w*), and how high in the boundary layer the measurement is (z/zi). Crucially, the correction is larger when convection is strong (hot, sunny, low-wind conditions) and smaller when wind-driven mixing dominates. It also predicts that sensible heat flux is proportionally more underestimated than latent heat flux (ET), which makes physical sense given how temperature and humidity are distributed vertically in the atmosphere.

- What this paper adds: This study takes that correction and actually applies it to real-world data from three European sites — a beech forest in Denmark (tall, windy, relatively good closure to begin with) and two grassland sites in the Bavarian pre-Alps (short vegetation, low measurement height, worse closure). The key question: does the correction work in practice?

- The results: For the Danish forest site, where the tower was tall enough that the correction could be applied directly, it nearly perfectly closed the energy balance — bringing the regression slope of turbulent fluxes versus available energy from 0.94 to 0.99. For the two grassland sites, the measurement height was too low for the correction to be applied directly (the LES didn't have fine enough resolution below about 20 m), so it had to be scaled using the observed daily energy balance ratio. Even so, the correction substantially improved things.

- The most compelling validation came from lysimeters — instruments that directly measure ET by weighing soil-filled tanks and tracking water loss. When the authors compared ET estimated from the corrected eddy-covariance data against the lysimeter measurements, the LES-based correction outperformed every other method tested. The typical error (RMSE) was cut roughly in half compared to leaving the data uncorrected, and the bias dropped to nearly zero (under 0.1 mm/day in absolute terms). In contrast, the most common previous method — assuming the energy imbalance should be split between sensible and latent heat in proportion to the Bowen ratio — consistently overestimated ET.