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

### Step 1: Implement Your ODE System
Edit the `ode` function in `main.cpp`:  
`auto ode(double t, const std::vector<double>& y, const Params& p) -> std::vector<double>`
//Where vector y is where you enter your initial conditions (position, velocity, acceleration) and Params p is where you enter your constants (mass, gravity, length, etc)  
Example: `y0 = {x,x_p}` (2 degrees of freedom, position and velocity)   
//Add the O.D.E. logic inside the function body  
 Example: $\ddot{x} + \frac{m}{L} x = 0$  
`return { std::vector<double> };` -> Returns a vector with the logic applied.   
// The returned vector must be the same size and in the same order as the initial state vector    

### Step 2: Adjust Parameters  
`Params p = {/* new parameters /};`  
`std::vector<double> y0 = {/ new initial conditions */};`  
`double h = 0.1;` // New time step in seconds  
`double tf = 5.0;` // New final time in seconds  


## Example: Mass-Pendulum System
The code includes an example implementation for the pendulum attached to a moving body:  
$\ddot{\theta} = \frac{(m_1 + m_2)g \sin\theta - m_1 L \dot{\theta}^2 \sin\theta \cos\theta}{(m_1 + m_2)L - m_1 L \cos^2\theta}$
$\ddot{x}_2 = \frac{m_1 L \dot{\theta}^2 \sin\theta - m_1 g \sin\theta \cos\theta}{(m_1 + m_2)L - m_1 L \cos^2\theta}$


### Default Parameters
`Params p = {1.0, 2.0, 1.0, 9.81};` // $m₁, m₂, L, g$  
`std::vector<double> y0 = {M_PI/6, 0.0, 0.0, 0.0};` // $\theta, \dot\{\theta}, x_2, \dot\{x_2}$  



