#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>

typedef struct particle
{
  double pos_x;
  double pos_y;
  double mass;
  double vel_x;
  double vel_y;
  double brightness;
}particle_t;

int getInput(const char* __restrict filename, particle_t* __restrict data, const int N);
int setOutput(const char* __restrict filename, particle_t* __restrict data, const int N);
void calculateNewPositions(particle_t* particles, const int N, const int nsteps, const double dt);

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

int main(int argc, char* argv[])
{
  double time1 = get_wall_seconds();
  if (argc != 6)
  {
    printf("Check your input..\n");
    return -1;
  }
  const int N = atoi(argv[1]);
  const char* filename = argv[2];
  const int nsteps = atoi(argv[3]);
  const double dt = atof(argv[4]);
  //const int graphics = atoi(argv[5]);
  const char* filename_out = "result.gal";
  particle_t* data = (particle_t*)malloc(N*sizeof(particle_t));

  int err = getInput(filename, data, N);
  if(err == -1)
  {
    return -1;
  }
  calculateNewPositions(data, N, nsteps, dt);

  err = setOutput(filename_out, data, N);
  if(err == -1)
  {
    return -1;
  }

  free(data);
  printf("galsim main took %7.3f wall seconds.\n", get_wall_seconds()-time1);
  return 0;
}

int getInput(const char* __restrict filename, particle_t* __restrict data, const int N)
{
  FILE *stream_in;
  stream_in = fopen(filename,"rb");
  double* arr = (double*)malloc(N*6*sizeof(double));

  if(stream_in == NULL)
  {
    printf("Error: unable to open file: %s\n", filename);
    return -1;
  }
  size_t input_size = N*6*sizeof(double);

  /* Read input to data array */
  size_t items_read = fread(arr, sizeof(char), input_size, stream_in);
  int j = 0;
  for (int i = 0; i < N*6; i += 6)
  {
    data[j].pos_x = arr[i];
    data[j].pos_y = arr[i+1];
    data[j].mass = arr[i+2];
    data[j].vel_x = arr[i+3];
    data[j].vel_y = arr[i+4];
    data[j].brightness = arr[i+5];
    j++;
  }

  if (items_read != input_size)
  {
    printf("Error reading the input file.\n");
    return -1;
  }

  fclose(stream_in);
  free(arr);
  return 1;
}

int setOutput(const char* __restrict filename, particle_t* __restrict data, const int N)
{
  FILE *stream_out;
  stream_out = fopen(filename,"wb");
  if(stream_out == NULL)
  {
    printf("Error: unable to open file: %s\n", filename);
    return -1;
  }
  fwrite(data, sizeof(particle_t), N, stream_out);
  fclose(stream_out);
  return 0;
}

void calculateNewPositions(particle_t* particles, const int N, const int nsteps, const double dt)
{
  double Fx, Fy,r, ax, ay;
  const double G = 100/(double)(N);
  const double epsilon = 0.001;
  particle_t *particles_new = (particle_t*)malloc(N*sizeof(particle_t));
  for (int n = 0; n < nsteps; n++)
  {
    for(int i = 0; i < N; i++)
    {
      Fx = 0;
      Fy = 0;
      for(int j = 0; j < i; j++)
      {
          r = sqrt((particles[i].pos_x - particles[j].pos_x)*(particles[i].pos_x - particles[j].pos_x) + (particles[i].pos_y - particles[j].pos_y)*(particles[i].pos_y - particles[j].pos_y));
          Fx += particles[j].mass*(particles[i].pos_x - particles[j].pos_x)/((r+epsilon)*(r+epsilon)*(r+epsilon));
          Fy += particles[j].mass*(particles[i].pos_y - particles[j].pos_y)/((r+epsilon)*(r+epsilon)*(r+epsilon));
      }
      for(int j = i+1; j < N; j++)
      {
          r = sqrt((particles[i].pos_x - particles[j].pos_x)*(particles[i].pos_x - particles[j].pos_x) + (particles[i].pos_y - particles[j].pos_y)*(particles[i].pos_y - particles[j].pos_y));
          Fx += particles[j].mass*(particles[i].pos_x - particles[j].pos_x)/((r+epsilon)*(r+epsilon)*(r+epsilon));
          Fy += particles[j].mass*(particles[i].pos_y - particles[j].pos_y)/((r+epsilon)*(r+epsilon)*(r+epsilon));
      }
      /* Forces */
      Fx *= -G*particles[i].mass;
      Fy *= -G*particles[i].mass;
      /* Acceleration */
      ax = Fx/(particles[i].mass);
      ay = Fy/(particles[i].mass);
      /* Updating velocities */
      particles_new[i].vel_x = particles[i].vel_x + dt*ax;
      particles_new[i].vel_y = particles[i].vel_y + dt*ay;
      /* Updating positions */
      particles_new[i].pos_x = particles[i].pos_x + dt* particles_new[i].vel_x;
      particles_new[i].pos_y = particles[i].pos_y + dt* particles_new[i].vel_y;
    }
    for (int i = 0; i < N; i++)
    {
      particles[i].pos_x = particles_new[i].pos_x;
      particles[i].pos_y = particles_new[i].pos_y;
      particles[i].vel_x = particles_new[i].vel_x;
      particles[i].vel_y = particles_new[i].vel_y;
    }
  }
  free(particles_new);
}
