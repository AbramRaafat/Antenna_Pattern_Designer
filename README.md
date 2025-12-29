# Antenna Pattern Designer

A collection of MATLAB-based tools for antenna design and analysis. This repository currently features the **Aperture Blockage App**, a tool designed to simulate and analyze the effects of physical blockages on rectangular antenna apertures.

The tools are built using MATLAB App Designer, providing an interactive GUI for visualization of geometry, radiation patterns (2D & 3D), and performance metrics.

---

##  App 1: Aperture Blockage Analyzer

This application simulates a rectangular aperture antenna with a rectangular blockage (e.g., a feed horn or sub-reflector strut). It calculates the degradation in gain, side lobe levels (SLL), and directivity using the **Principle of Superposition**.

### Key Features
* **Parametric Geometry:** Fully adjustable aperture and blockage dimensions.
* **Scanning Capabilities:** Simulates beam scanning effects (Scanning H-Plane).
* **Multi-View Analysis:** Dedicated tabs for H-Plane, E-Plane, Polar, and 3D views.
* **Performance Metrics:** Automatic calculation of **Gain Degradation**, **Max SLL**, **Directivity**, and **Effective Area**.

---

##  Walkthrough & Tutorial

### 1. Defining the Geometry
When you launch the app, you can define the main aperture dimensions (`a`, `b`) and the blockage dimensions (`cx`, `cy`). In the **Geometry Tab**, the blue rectangle represents the clear aperture, and the red rectangle represents the blockage.

![Geometry Center](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/A3.png)
*Figure 1: A standard setup with a central blockage ($c_x=1, c_y=2$) inside a $3\lambda \times 4\lambda$ aperture.*

### 2. Adjusting Blockage Configuration
The app allows you to explore different blockage scenarios, such as struts or offset feeds.
* **Geometry Tab:** Visualizes the physical layout.
* **Inputs:** Use `Offset X` and `Offset Y` to simulate asymmetric blockage.

![Geometry Offset](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/A4.png)
*Figure 2: An offset blockage configuration, useful for analyzing feed offset effects.*

### 3. analyzing 2D Cuts (H-Plane, E-Plane, Polar)
Before viewing the full 3D model, the app provides detailed 2D cuts to analyze specific beam characteristics:

* **H-Plane Tab (Red Curve):**
    * Displays the pattern along the scanning plane ($u = \sin(\theta)$).
    * Includes a vertical dashed line indicating the current **Scan Angle**.
    * Useful for observing grating lobes or asymmetry caused by scanning.

* **E-Plane Tab (Blue Curve):**
    * Displays the pattern along the transverse (orthogonal) plane.
    * Shows how the blockage affects the beam width in the non-scanning dimension.

* **Polar Tab:**
    * Overlays both the H-Plane (Red) and E-Plane (Blue) cuts on a polar coordinate system.
    * Provides a direct visual comparison of the **Beamwidth** and **Side Lobe Levels** between the two planes.

### 4. 3D Radiation Pattern
The **3D Pattern Tab** combines all data into an interactive surface plot. This allows you to see the "total" effect of the blockage, including diagonal side lobes that might not appear in the principal cuts.

* **Red Line:** Traces the H-Plane cut.
* **Blue Line:** Traces the E-Plane cut.

![3D Pattern](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/A1.png)
*Figure 3: The final 3D representation showing the beam magnitude. The results panel on the left updates instantly to show a Gain Degradation of -3.52 dB and an Effective Area of 8.00 sq.lam.*

---

## Theoretical Background

The app uses the **Field Equivalence Principle** and **Pattern Multiplication**. The total field is calculated by subtracting the field of the blockage from the field of the ideal aperture:

$$E_{total} = E_{aperture} - E_{blockage}$$

Where the far-field pattern of a rectangular aperture is approximated using Sinc functions:

$$E(\theta) \propto Area \cdot \text{sinc}(k \cdot \text{width} \cdot u)$$

* **Scanning:** The app introduces a phase shift $\gamma = \beta \sin(\theta_0)$ to steer the beam.
* **Blockage Phase:** If the blockage is offset by $(d_x, d_y)$, a phase term $e^{j2\pi(d_x u + d_y v)}$ is applied to the subtracted blockage pattern.

---
---

##  App 2: Linear Dipole Analyzer

This tool provides a deep dive into the physics of linear wire antennas. It visualizes the relationship between the physical length of the dipole, the current distribution along the wire, and the resulting far-field radiation pattern.

### Key Features
* **Variable Length Analysis:** Simulate dipoles from short fractions of a wavelength up to long wires ($10\lambda$).
* **Current Distribution Models:** Compare **Sinusoidal** (Real), **Triangular** (Short approximation), and **Uniform** (Hertzian ideal) distributions.
* **Impedance Calculation:** Estimates Radiation Resistance ($R_{rad}$) based on the current integrals.
* **Visual Validation:** Alerts users when specific approximations (like the Short Dipole model) become invalid for the selected length.

###  Walkthrough

#### 1. Configuration & Modeling
In the control panel, you can set the dipole length in terms of wavelengths ($L/\lambda$). You can also select the mathematical model for the current distribution.

![Dipole Configuration](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/D3.png)
*Figure 5: Setting the length to $0.5\lambda$ (Half-Wave Dipole) and selecting the General Sinusoidal model.*

#### 2. Visualizing Current Distribution
The **Current Distribution Tab** shows how the current magnitude varies along the z-axis.
* **Blue Line:** The normalized current magnitude $|I(z)|$.
* **Physics Note:** Notice how the current naturally goes to zero at the ends of the wire (boundary conditions). For very short dipoles, this looks triangular.

![Current Distribution](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/D2.png)
*Figure 6: A triangular current distribution, which is a valid approximation for short dipoles ($L \ll \lambda$).*

#### 3. Radiation Pattern Analysis
The app calculates the pattern by integrating the current distribution.
* **2D Polar Plot:** Shows the E-Plane cut (the classic "figure-eight" for a half-wave dipole).
* **3D Plot:** Visualizes the full "donut" shape characteristic of omnidirectional wire antennas.
* **Metrics:** The dashboard displays calculated **Directivity** (e.g., 1.76 dBi for a short dipole, 2.15 dBi for $\lambda/2$) and **Half-Power Beamwidth (HPBW)**.

![Dipole Pattern](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/D1.jpg)
*Figure 7: The radiation pattern view showing the omnidirectional nature of the dipole in 3D and the elevation cut in 2D.*

###  Theoretical Background
The far-field pattern $F(\theta)$ is derived from the finite integration of the current $I(z)$:

$$F(\theta) = \frac{\cos(kH \cos \theta) - \cos(kH)}{\sin \theta}$$

Where $H = L/2$. The app numerically integrates the Poynting vector to find the total radiated power ($P_{rad}$) and Radiation Resistance:

$$R_{rad} = \frac{2 P_{rad}}{|I_{max}|^2}$$

---

## 📡 App 3: Parabolic Reflector Designer

A high-frequency design tool for reflector antennas. This app calculates the efficiency breakdown (Spillover, Taper, Blockage) and simulates the focusing properties of parabolic dishes using Ray Tracing and aperture integration.

### Key Features
* **Dual Configurations:** Supports both **Front-Fed** and **Cassegrain** (dual-reflector) designs.
* **Feed Modeling:** Select between an **Ideal Horn** (cosine power $q$) or a **Dipole Feed**.
* **Efficiency Budget:** Automatically calculates Total Aperture Efficiency ($\epsilon_{ap}$) based on spillover, taper, and blockage.
* **Ray Tracing:** Visualizes the focal point and sub-reflector ray paths.

### 📖 Walkthrough

#### 1. System Design & Ray Tracing
The **Ray Tracing Tab** provides a geometric verification of your design. It visualizes the primary dish, the focal point, and the feed blockage.
* **Inputs:** Define Frequency, Dish Diameter, and $f/D$ ratio.
* **Interactive Sync:** Adjusting the $f/D$ ratio automatically updates the subtended angle ($\theta_0$), and vice versa.

![Reflector Geometry](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/P1.png)
*Figure 8: A Front-Fed configuration showing rays focusing at the feed point. The red lines represent the blockage shadow.*

#### 2. Feed Selection & Polarization
You can model different feed sources to see how they illuminate the dish.
* **Controls:** Select "Dipole (Eq. 16)" or "Horn" and adjust the horn's cosine power factor ($q$).
* **Polarization Tab:** Visualizes the field vector alignment on the aperture. Note that a simple dipole feed creates cross-polarization (curved field lines) compared to an ideal Huygens source.

| Feed Configuration | Polarization Analysis |
| :---: | :---: |
| ![Feed Config](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/P3.png) | ![Polarization](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/P5.png) |
| *Figure 9: Selecting a Dipole Feed.* | *Figure 10: The resulting Field lines on the aperture.* |

#### 3. Performance Metrics & Efficiency
Reflector design is a trade-off between **Spillover** (energy missing the dish) and **Taper** (uneven illumination). The app calculates these instantly:

![Design Metrics](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/P4.png)
*Figure 11: The metrics panel showing a Directivity of 43.34 dBi. The efficiency breakdown shows how the 60.8% total efficiency is derived from Spillover, Taper, and Blockage losses.*

#### 4. 3D Far-Field Pattern
Finally, the app performs a numerical integration of the aperture field (using Bessel functions for the circular aperture) to generate the high-gain pencil beam pattern.

![3D Pencil Beam](https://github.com/AbramRaafat/Antenna_Pattern_Designer/blob/main/imgs/P6.png)
*Figure 12: 3D visualization of the main beam and the first few side lobes.*

###  Theoretical Background
The Directivity $D$ is calculated based on the effective area $A_{eff}$:

$$D = \frac{4\pi}{\lambda^2} A_{phys} \cdot \epsilon_{ap}$$

Where the total efficiency $\epsilon_{ap}$ is the product of sub-efficiencies:

$$\epsilon_{ap} = \epsilon_{spillover} \cdot \epsilon_{taper} \cdot \epsilon_{phase} \cdot \epsilon_{blockage}$$

The app calculates blockage efficiency by subtracting the area of the feed/sub-reflector:
$$\epsilon_{block} = (1 - (\frac{D_{block}}{D_{main}})^2)^2$$

