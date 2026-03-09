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

- An earlier study (**De Roo et al., 2018**) used large-eddy simulations — very detailed computer models of the lower atmosphere — to figure out how large these missed dispersive fluxes are under different atmospheric conditions. The result was a set of correction equations that take in two key variables: how strong the surface friction is relative to the buoyancy-driven convection (expressed as u*/w*), and how high in the boundary layer the measurement is (z/zi). Crucially, the correction is larger when convection is strong (hot, sunny, low-wind conditions) and smaller when wind-driven mixing dominates. It also predicts that sensible heat flux is proportionally more underestimated than latent heat flux (ET), which makes physical sense given how temperature and humidity are distributed vertically in the atmosphere.

- What this paper adds: This study takes that correction and actually applies it to real-world data from three European sites — a beech forest in Denmark (tall, windy, relatively good closure to begin with) and two grassland sites in the Bavarian pre-Alps (short vegetation, low measurement height, worse closure). The key question: does the correction work in practice?

- The results: For the Danish forest site, where the tower was tall enough that the correction could be applied directly, it nearly perfectly closed the energy balance — bringing the regression slope of turbulent fluxes versus available energy from 0.94 to 0.99. For the two grassland sites, the measurement height was too low for the correction to be applied directly (the LES didn't have fine enough resolution below about 20 m), so it had to be scaled using the observed daily energy balance ratio. Even so, the correction substantially improved things.

- The most compelling validation came from lysimeters — instruments that directly measure ET by weighing soil-filled tanks and tracking water loss. When the authors compared ET estimated from the corrected eddy-covariance data against the lysimeter measurements, the LES-based correction outperformed every other method tested. The typical error (RMSE) was cut roughly in half compared to leaving the data uncorrected, and the bias dropped to nearly zero (under 0.1 mm/day in absolute terms). In contrast, the most common previous method — assuming the energy imbalance should be split between sensible and latent heat in proportion to the Bowen ratio — consistently overestimated ET.

**Stoughton et al. 2002**

- What's this study about? Researchers spent two weeks in summer 1996 measuring how much water evaporates from a pecan orchard in the New Mexico desert near Las Cruces. They wanted to understand how the dry surrounding desert affects water loss from the irrigated trees, and how air moves through the orchard.

- Why does it matter? Irrigated farms in deserts use enormous amounts of water. Accurately predicting how much water an orchard actually consumes is critical for water management in arid regions. The problem is that standard calculation methods weren't designed for orchards sitting in a desert surrounded by hot, dry air.

- What makes this orchard unusual?
The pecan trees only covered about half the ground — the rest was bare, sun-baked soil between the tree rows. This open structure made the orchard behave very differently from a dense forest or a solid field of crops.

Key findings:
- The orchard used significantly more water than standard formulas predicted. The most common method (Penman-Monteith) underestimated actual water use by 82%.
- A more sophisticated method that accounts for the "oasis effect" — hot desert air blowing over the irrigated orchard and adding extra energy — was much closer, underestimating by only 11%.
- The hot bare soil between the tree rows acted like a furnace, pushing warm air upward through the open canopy spaces. This made the orchard behave like a "weak oasis" — somewhat cooler and wetter than the surrounding desert, but not as extreme as a fully irrigated, dense crop would be.
- Air moved very freely up and down through the orchard, unlike in dense forests where the canopy blocks vertical airflow.

- Bottom line: Standard evaporation formulas significantly undercount water use in open desert orchards because they don't adequately capture the influence of surrounding hot, dry air. Water managers should use more advanced approaches when planning irrigation for orchards in arid climates.

**Alexander et al (2022)** Simulating land-atmosphere coupling in the Central Valley, California: Investigating soil moisture impacts on boundary layer properties

- Researchers used a computer weather model to simulate summer conditions in California's Central Valley (CV) — one of the most intensively irrigated agricultural regions on Earth — to figure out why these simulations often produce inaccurate results. Specifically, they wanted to understand how the moisture in the soil affects the lowest layer of the atmosphere, and which model settings best capture what's actually happening.

- The Central Valley produces over $20 billion worth of crops annually and relies almost entirely on irrigation since there's essentially no rain in summer. Accurate atmospheric models are essential for forecasting air quality, managing water use, predicting wildfire risk, and planning under climate change. Yet models consistently struggle to simulate conditions there correctly.

- They ran nine different versions of a standard weather model (WRF), each using different combinations of two key components: a "land surface model" (LSM) — which controls how the model handles soil moisture, evaporation, and plant water use — and a "planetary boundary layer scheme" (PBL) — which controls how the model mixes air near the ground. They compared outputs against a rich dataset of ground stations, aircraft flights, and weather instruments.

Key findings:
- The choice of land surface model mattered far more than the choice of boundary layer scheme. In other words, getting the soil and vegetation physics right is more important than getting the turbulence mixing right.
- All models struggled to capture the actual moisture in the system — they simulated air that was too dry and soil that was too dry, likely because the model starting conditions don't adequately represent how much irrigation is actually happening.
- No single model configuration performed best across all variables. One model (PX-NO) best captured latent heat and humidity, while others (Noah, Noah-MP) better reproduced temperatures, but for potentially the wrong reasons.
- A "smart" soil moisture correction scheme (PX with nudging) actually made things worse rather than better, likely because the reference data it was nudging toward was also too hot and dry for irrigated cropland.
- Potential evapotranspiration (ETo — the amount of water crops could theoretically use) was well-simulated across all models, but this was largely because ETo depends mostly on solar radiation, which models handle well. It doesn't mean the models correctly capture actual water use.

Bottom line: To improve weather and climate simulations of the Central Valley, the most urgent need is better soil moisture data — particularly from irrigated farmland — to fix the starting conditions models use. The models themselves also need updated physics to realistically represent irrigation. Without these improvements, the models are systematically too hot, too dry, and produce atmospheric layers that are unrealistically deep.


Sutton’s solution → problem

