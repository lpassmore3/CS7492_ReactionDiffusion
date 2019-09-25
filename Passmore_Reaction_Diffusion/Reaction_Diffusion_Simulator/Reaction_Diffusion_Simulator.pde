// CS 7492 - Simulation of Biology
// Project 3 - Reaction-Diffusion Simulator
// Author: Austin Passmore

boolean isSimulating = true;    // Controls if the simulation is paused or not
boolean isDisplayingU = true;    // Controls if u or v concentrations are drawn.
boolean isReacting = true;    // Controls if the reaction step occurs.
boolean isSpatiallyVarying = false;    // Controls whether k and f parameters vary spatially.

int frame = 1;    // Frame counter
int frameMax = 10;    // Ensures the screen is drawn every max frames.

int cellSize = 4;    // each grid space is 4x4 pixels
int gridSizeX = 200;  // the grid will be 200x200 cells
int gridSizeY = 200;

// grid[i][j], i == y and j == x
float[][] u = new float[gridSizeY][gridSizeX];    // Holds u concentrations
float[][] v = new float[gridSizeY][gridSizeX];    // Holds v concentrations

float dt = 2.5;        // Timestep
float ru = 0.082;    // Diffusion rate for u
float rv = 0.041;    // Diffusion rate for v
float k = 0.0625;
float f = 0.035;

void setup() {
  size(800, 800);
  background(0, 0, 0);
  noStroke();
  
  initConcentrations();    // Initialize the u and v concentrations specifically.
}

void draw() {
  
  if (isSimulating) {    // Checks if the simulation is paused or not
    updateConcentrations();    // Calls the diffusion and reaction steps  
  }
  
  // Ensures the screen is drawn every 10 frames.
  if (frame == frameMax) {
    drawCells();
    frame = 1;
  } else {
    frame++;
  }
}

// Controls keyboard inputs
void keyPressed() {
  if (key == 'u') {    // Display u concentrations
    isDisplayingU = true;
    println("Displaying u concentrations.");
  } else if (key == 'v') {    // Display v concentrations
    isDisplayingU = false;
    println("Displaying v concentrations.");
  } else if (key == 'i') {    // Reset screen to initial conditions.
    initConcentrations();
  } else if (key == 'd') {    // Toggle between diffusion alone and reaction-diffusion
    if (isReacting) {
      isReacting = false;
      println("Diffusion only mode");
    } else {
      isReacting = true;
      println("Reaction-Diffusion mode");
    }
  } else if (key == ' ') {    // Toggle between pausing and playing the simulation
    if (isSimulating) {
      isSimulating = false;
      println("Simulation paused.");
    } else {
      isSimulating = true;
      println("Resuming simulation.");
    }
  } else if (key == 'p') {    // Toggle between constant and spatially-varying parameters k and f
    if (isSpatiallyVarying) {
      isSpatiallyVarying = false;
      println("Constant parameter mode.");
    } else {
      isSpatiallyVarying = true;
      println("Spatially-varying parameter mode.");
    }
  } else if (key == '1') {    // Set parameters for spots
    k = 0.0625;
    f = 0.035;
    println("Parameters for spots set.");
  } else if (key == '2') {    // Set parameters for stripes
    k = 0.06;
    f = 0.035;
    println("Parameters for stripes set.");
  } else if (key == '3') {    // Set parameters for spiral waves
    k = 0.0475;
    f = 0.0118;
    println("Parameters for spiral waves set.");
  } else if (key == '4') {    // Set parameters for squiggle lines
    k = 0.055;
    f = 0.0275;
    println("Parameters for squiggle lines set.");
  }
}

// Prints the concentrations of cell that has been clicked on.
void mousePressed() {
  int iIndex = mouseY / cellSize;
  int jIndex = mouseX / cellSize;
  float uConcentration = u[iIndex][jIndex];
  float vConcentration = v[iIndex][jIndex];
  println("u =", uConcentration, ", v =", vConcentration);
}

// Colors each cell on a grayscale depending on min and max concentration.
// Display concentrations for u or v is determined by a flag with keyboard control.
void drawCells() {
  float min = 0;
  float max = 0;
  for (int i = 0; i < gridSizeY; i++) {
    for (int j = 0; j < gridSizeX; j++) {
      if (i == 0 && j == 0) {
        if (isDisplayingU) {
          min = u[i][j];
          max = u[i][j];
        } else {
          min = v[i][j];
          max = v[i][j];
        }
      } else {
        if (isDisplayingU) {
          if (u[i][j] < min) {
            min = u[i][j];
          }
          if (u[i][j] > max) {
            max = u[i][j];
          }
        } else {
          if (v[i][j] < min) {
            min = v[i][j];
          }
          if (v[i][j] > max) {
            max = v[i][j];
          }
        }
      }
    }
  }
  for (int i = 0; i < gridSizeY; i++) {
    for (int j = 0; j < gridSizeX; j++) {
      float concentration = 0;
      if (isDisplayingU) {
        concentration = u[i][j];
      } else {
        concentration = v[i][j];
      }
      float grayscale = ((concentration - min) / (max - min)) * 255;
      fill(grayscale, grayscale, grayscale);
      rect(j * cellSize, i * cellSize, cellSize, cellSize);
    }
  }
}

// Initializes all cells to have u = 1.0, v = 0.0
// and have a 10x10 square to have u =0.5, v= 0.25
public void initConcentrations() {
  for (int i = 0; i < gridSizeY; i++) {
    for (int j = 0; j < gridSizeX; j++) {
      // Create the 10x10 block in the center of the screen
      if (i >= (gridSizeY / 2) - 5 && i <= (gridSizeY / 2) + 5 && j >= (gridSizeX / 2) - 5 && j <= (gridSizeX / 2) + 5) {
        u[i][j] = 0.5;
        v[i][j] = 0.25;
      } else {
        u[i][j] = 1.0;
        v[i][j] = 0.0;
      }
    }
  }
}

// Updates the concentration value by performing the diffusion and reaction steps.
public void updateConcentrations() {
  diffusionStep(u, dt, ru);
  diffusionStep(v, dt, rv);
  if (isReacting) {
    reactionStep(u, v, dt);
  }
}

// Performs diffusion of a single element (u or v) with a hard boundary by calculating the Laplacian.
public void diffusionStep(float[][] x, float dt, float rx) {
  float[][] x0 = x;
  for (int i = 0; i < gridSizeY; i++) {
    for (int j = 0; j < gridSizeX; j++) {
      float laplacian;
      if (i == 0 && j == 0) {    // Top-left corner
        laplacian = x0[i][j + 1] + x0[i + 1][j] - (2 * x0[i][j]);
      } else if (i == 0 && j == gridSizeX - 1) {    // Top-right corner
        laplacian = x0[i][j - 1] + x0[i + 1][j] - (2 * x0[i][j]);
      } else if (i == gridSizeY - 1 && j == 0) {    // Bottom-left corner
        laplacian = x0[i - 1][j] + x0[i][j + 1] - (2 * x0[i][j]);
      } else if (i == gridSizeY - 1 && j == gridSizeX - 1) {    // Bottom-right corner
        laplacian = x0[i - 1][j] + x0[i][j - 1] - (2 * x0[i][j]);
      } else if (i == 0 && j > 0 && j < gridSizeX - 1) {    // Top boundary
        laplacian = x0[i][j - 1] + x0[i][j + 1] + x0[i + 1][j] - (3 * x0[i][j]); 
      } else if (i == gridSizeY - 1 && j > 0 && j < gridSizeX - 1) {    // Bottom boundary
        laplacian = x0[i][j - 1] + x0[i][j + 1] + x0[i - 1][j] - (3 * x0[i][j]); 
      } else if (i > 0 && i < gridSizeY - 1 && j == 0) {    // Left boundary
        laplacian = x0[i - 1][j] + x0[i + 1][j] + x0[i][j + 1] - (3 * x0[i][j]); 
      } else if (i > 0 && i < gridSizeY - 1 && j == gridSizeX - 1) {    // Right boundary
        laplacian = x0[i - 1][j] + x0[i + 1][j] + x0[i][j - 1] - (3 * x0[i][j]); 
      } else {
        laplacian = x0[i - 1][j] + x0[i + 1][j] + x0[i][j - 1] + x0[i][j + 1] - (4 * x0[i][j]); 
      }
      float dx = rx * dt * laplacian;
      x[i][j] += dx;
    }
  }
}

// Performs the reaction between both elements u and v.
// u += -uv^2 + f(1 - u)
// v += uv^2 - (f + k)v
public void reactionStep(float[][] u, float[][] v, float dt) {
  for (int i = 0; i < gridSizeY; i++) {
    for (int j = 0; j < gridSizeX; j++) {
      if (isSpatiallyVarying) {
        k = 0.03 + (0.04 * j / (gridSizeX - 1));
        f = 0.08 - (0.08 * i / (gridSizeY - 1));
      }
      float uValue = u[i][j];
      float vValue = v[i][j];
      float uv2 = uValue * vValue * vValue;
      u[i][j] += ((-1 * uv2) + (f * (1 - uValue))) * dt;
      v[i][j] += (uv2 - ((f + k) * vValue)) * dt;
    }
  }
}
