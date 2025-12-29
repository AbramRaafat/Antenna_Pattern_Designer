# Antenna Pattern Designer

A collection of MATLAB-based tools for antenna design and analysis. This repository currently features the **Aperture Blockage App**, a tool designed to simulate and analyze the effects of physical blockages on rectangular antenna apertures.

The tools are built using MATLAB App Designer, providing an interactive GUI for visualization of geometry, radiation patterns (2D & 3D), and performance metrics.

---

## 📡  Aperture Blockage Analyzer

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


