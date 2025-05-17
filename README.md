# Generic ODE Solver with RK4

A generic numerical solver for systems of Ordinary Differential Equations (ODEs) using the 4th order Runge-Kutta (RK4) method. Includes an example of a coupled mass-pendulum system.

## Main Features
✅ **Generic ODE solving**  
Adaptable to any system of first-order differential equations by editing the `ode` function.

⏱️ **Benchmark Mode**  
Accurately measures execution time for performance analysis.

📊 **CSV Output**  
Automatically generates a file with the simulation results.

🔧 **Modular Structure**  
Easily implement new physical or mathematical systems.

---

## Usage

### 1. Building
`git clone https://github.com/JeeerryZ/ode-solver.git`  
`cd ode-solver`  
`cmake -B build -S .`  
`cmake --build build --config Release`

### 2. Basic Execution
`./ode-solver` # Generates results.csv with simulation data.

### 3. Benchmark algorithm
`./ode-solver --benchmark` # Shows time spent calculating your ode.

### 4. Custom Output File
`./ode-solver --output myfile` # At the moment we only support .csv files so you dont need to specify it.


## Adapting to Other Systems

## Step 1: Implement Your ODE System

Edit the `ode` function in `main.cpp` following this template:
```cpp
auto ode(double t, const std::vector<double>& y, const Params& p) -> std::vector<double> {  
    // Your ODE logic here  
}
```

Parameter Structure

| Parameter | Description                              | Example Usage      |
|-----------|------------------------------------------|--------------------|
| y         | State vector (current values)            | y = position, y = velocity |
| p         | System constants (Params struct)         | p.m1, p.L, etc.    |

Example: Simple Harmonic Oscillator

Equation:
$$
\ddot{x} + \frac{m}{L} x = 0
$$

Implementation:
```cpp
auto ode(double t, const std::vector<double>& y, const Params& p) -> std::vector<double> {
    double x = y;      // Position
    double x_p = y[1];    // Velocity

    double x_pp = -p.m / p.L * x;  // Acceleration

    return {x_p, x_pp};  // Returns [dx/dt, d²x/dt²]
}
```
Essential Rules

1. Dimensional consistency:
   The returned vector must have the same size and order as the input state vector.  

   Example:  
   // Input state: [θ, θ']  
   // Correct return: [θ', θ'']  
   return {theta_p, theta_pp};

### Step 2: Adjust Parameters  
```cpp
Params p{} //Add your constants
std::vector<double> y0{} // Set your initial conditions
double h = 0.1; // Set your time step
double tf = 5.0; // Set the total time of the simulation
```


## Example: Mass-Pendulum System
The code includes an example implementation for the pendulum attached to a moving body:  
$\ddot{\theta} = \frac{(m_1 + m_2)g \sin\theta - m_1 L \dot{\theta}^2 \sin\theta \cos\theta}{(m_1 + m_2)L - m_1 L \cos^2\theta}$
$\ddot{x}_2 = \frac{m_1 L \dot{\theta}^2 \sin\theta - m_1 g \sin\theta \cos\theta}{(m_1 + m_2)L - m_1 L \cos^2\theta}$


### Default Parameters
```cpp
Params p = {1.0, 2.0, 1.0, 9.81}; // $m₁, m₂, L, g$  
std::vector<double> y0 = {M_PI/6, 0.0, 0.0, 0.0}; // $\theta, \dot\{\theta}, x_2, \dot\{x_2}$  
```


