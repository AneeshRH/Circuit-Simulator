---

# **ARHSpice – Circuit Simulator**

![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)
![Status](https://img.shields.io/badge/Build-Passing-brightgreen.svg)
![License](https://img.shields.io/badge/License-MIT-green.svg)
![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20Mac-lightgrey)

ARHSpice is a lightweight SPICE-style circuit simulator with a drag-and-drop GUI and a symbolic Laplace-domain Modified Nodal Analysis (MNA) solver.
It allows you to **draw**, **simulate**, and **visualize** analog circuits entirely in Python.

---

## ⭐ Features

### 🎛️ Graphical Circuit Editor

* Drag-and-drop components
* Resistor, Capacitor, Inductor, Voltage Source, Current Source
* Ground placement
* Rotate, move, delete
* Wiring with automatic node merging
* Clean UI built using Tkinter

### 🧮 Symbolic Laplace-Domain Solver

* Full MNA formulation in the Laplace domain
* Supports **R, L, C, V, I**
* Automatic inverse Laplace transform
* Time-domain voltage & current waveforms
* Designed using pure Python + SymPy

### 📊 Plotting & Visualization

* Voltage and current graph tabs
* Multi-subplot plotting
* Zoom/pan using Matplotlib toolbar
* Adjustable simulation time range
* Checkboxes to select components to plot

---

## 📥 Installation

### 1. Install Dependencies

```bash
pip install sympy numpy matplotlib pillow networkx
```

### 2. Run the Application

```bash
python trial1.py
```

---

## 🧠 How It Works

### 1. Create the Circuit

Use the GUI to place components and wire them.

### 2. Generate Netlist

Click **SIMULATE** → creates `output.txt`

Example netlist:

```
R1 1 2 100
C1 2 0 10u
V1 1 0 5
```

### 3. Symbolic Solver (graphz.py)

Each component is converted to its Laplace admittance:

| Component | Laplace Model             |
| --------- | ------------------------- |
| R         | 1/R                       |
| C         | sC                        |
| L         | 1/(sL)                    |
| Vsrc      | Adds extra MNA row/column |
| Isrc      | RHS = I/s                 |

Then the solver performs:

```
[G(s)] · X(s) = I(s)
```

And returns the time-domain result using:

```
x(t) = L⁻¹{X(s)}
```

### 4. Plotting

Choose components → Set time range → Click **Update Plot**.

---

## 📁 Project Structure

```
ARHSpice/
│
├── trial1.py          # Main GUI (schematic editor + plotting)
├── graphz.py          # Symbolic Laplace-domain MNA solver
├── output.txt         # Generated netlist
│
├── resistor.png
├── capacitor.png
├── inductor.png
├── voltage_source.png
├── current_source.png
└── ground.png         # Component assets
```

---

## ⚠️ Limitations

* Only linear components (RLC + ideal sources)
* Voltage and current sources act as **step inputs**
* Large circuits may be slow due to symbolic inverse Laplace
* No AC/sinusoidal sources yet

---

## 📜 License

Released under the **MIT License**.
You are free to modify and distribute.

---

## 💡 Roadmap

* [ ] Add frequency-domain (AC) analysis
* [ ] Numerical transient solver (RK45)
* [ ] Support for diodes, BJTs, MOSFETs
* [ ] Import/export standard SPICE netlists
* [ ] Dark theme GUI
* [ ] Export plot images

---

## 🤝 Contributing

Pull requests and suggestions are welcome.
For major changes, please open an issue first.

---

