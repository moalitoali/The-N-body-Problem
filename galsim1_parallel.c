#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <sys/time.h>
#include <pthread.h>

typedef struct particle
{
  double pos_x;
  double pos_y;
  double mass;
  double vel_x;
  double vel_y;
  double brightness;
}particle_t;

typedef struct dataForThread
{
  int start;
  int stop;
  int N;
  double dt;
  pthread_t thread;
}dataForThread_t;

/* Global arrays */
particle_t* particles;
particle_t* particles_new;
int rest;

int getInput(const char* __restrict filename, const int N);
int setOutput(const char* __restrict filename, const int N);
void calculateNewPositions(const int N, const int nsteps, const double dt, int num_of_threads);
void* thread_func(void* arg);

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

int main(int argc, char* argv[])
{
  double time1 = get_wall_seconds();
  if (argc != 8)
  {
    printf("Expected input: ./galsim_A3 N filename nsteps delta_t theta_max graphics n_threads\n");
    return -1;
  }
  const int N = atoi(argv[1]);
  const char* filename = argv[2];
  const int nsteps = atoi(argv[3]);
  const double dt = atof(argv[4]);
  //const double theta_max = atof(argv[5]);
  //const int graphics = atoi(argv[6]);
  const int num_of_threads = atoi(argv[7]);

  rest = N%num_of_threads;

  const char* filename_out = "result.gal";
  particles = (particle_t*)malloc(N*sizeof(particle_t));

  int err = getInput(filename, N);
  if(err == -1)
  {
    return -1;
  }

  calculateNewPositions(N, nsteps, dt, num_of_threads);

  err = setOutput(filename_out, N);
  if(err == -1)
  {
    return -1;
  }

  free(particles);
  printf("galsim main took %7.3f wall seconds.\n", get_wall_seconds()-time1);
  return 0;
}

int getInput(const char* __restrict filename, const int N)
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
    particles[j].pos_x = arr[i];
    particles[j].pos_y = arr[i+1];
    particles[j].mass = arr[i+2];
    particles[j].vel_x = arr[i+3];
    particles[j].vel_y = arr[i+4];
    particles[j].brightness = arr[i+5];
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

int setOutput(const char* __restrict filename, const int N)
{
  FILE *stream_out;
  stream_out = fopen(filename,"wb");
  if(stream_out == NULL)
  {
    printf("Error: unable to open file: %s\n", filename);
    return -1;
  }
  fwrite(particles, sizeof(particle_t), N, stream_out);
  fclose(stream_out);
  return 0;
}

void calculateNewPositions(const int N, const int nsteps, const double dt, int num_of_threads)
{
  dataForThread_t threads[num_of_threads];
  particles_new = (particle_t*)malloc(N*sizeof(particle_t));
  for (int n = 0; n < nsteps; n++)
  {
    int work_size = N/num_of_threads; // number of particles for every thread
    /* Create threads */
    for(int i = 0; i < num_of_threads-1; i++){
      threads[i].start = (i+1) * work_size - work_size;
      threads[i].stop = threads[i].start + work_size;
      threads[i].N = N;
      threads[i].dt = dt;

      pthread_create(&(threads[i].thread), NULL, thread_func, &threads[i]);
    }
    /* Create last thread (with possibly different work_size) */
    threads[num_of_threads-1].start = num_of_threads * work_size - work_size;
    threads[num_of_threads-1].stop = threads[num_of_threads-1].start + work_size + rest;
    threads[num_of_threads-1].N = N;
    threads[num_of_threads-1].dt = dt;
    pthread_create(&(threads[num_of_threads-1].thread), NULL, thread_func, &threads[num_of_threads-1]);

    /* Join threads */
    for(int i = 0; i < num_of_threads; i++){
      pthread_join(threads[i].thread, NULL);
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

void* thread_func(void* arg){
  dataForThread_t* info = (dataForThread_t *) arg;

  double Fx, Fy, r, ax, ay;
  const double G = 100/(double)(info->N);
  const double epsilon = 0.001;
  /* Compute new position of particle */
  for(int i = info->start; i<info->stop; i++){ // for specific particles
    Fx = 0;
    Fy = 0;
    for(int j = 0; j<i; j++){ // compare with every other particle
      /* Distance between particles */
      r = sqrt((particles[i].pos_x - particles[j].pos_x)*(particles[i].pos_x - particles[j].pos_x) + (particles[i].pos_y - particles[j].pos_y)*(particles[i].pos_y - particles[j].pos_y));

      /* Forces*/
      Fx += particles[j].mass*(particles[i].pos_x - particles[j].pos_x)/((r+epsilon)*(r+epsilon)*(r+epsilon));
      Fy += particles[j].mass*(particles[i].pos_y - particles[j].pos_y)/((r+epsilon)*(r+epsilon)*(r+epsilon));
    }
    for(int j = i+1; j < info->N; j++){
      /* Distance between particles */
      r = sqrt((particles[i].pos_x - particles[j].pos_x)*(particles[i].pos_x - particles[j].pos_x) + (particles[i].pos_y - particles[j].pos_y)*(particles[i].pos_y - particles[j].pos_y));

      /* Forces*/
      Fx += particles[j].mass*(particles[i].pos_x - particles[j].pos_x)/((r+epsilon)*(r+epsilon)*(r+epsilon));
      Fy += particles[j].mass*(particles[i].pos_y - particles[j].pos_y)/((r+epsilon)*(r+epsilon)*(r+epsilon));
    }
    /* Forces */
    Fx *= -G * particles[i].mass;
    Fy *= -G * particles[i].mass;
    /* Acceleration */
    ax = Fx/(particles[i].mass);
    ay = Fy/(particles[i].mass);
    /* Updating velocities */
    particles_new[i].vel_x = particles[i].vel_x + info->dt * ax;
    particles_new[i].vel_y = particles[i].vel_y + info->dt * ay;
    /* Updating positions */
    particles_new[i].pos_x = particles[i].pos_x + (double)info->dt * particles_new[i].vel_x;
    particles_new[i].pos_y = particles[i].pos_y + (double)info->dt * particles_new[i].vel_y;
  }
  return NULL;
}