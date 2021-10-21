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
}__attribute__((__packed__)) particle_t;

typedef struct quadrant
{
  /* Position of specific quadrant in the whole space */
  double origin_x;
  double origin_y;

  double width;
  double height;
} __attribute__((__packed__)) quadrant_t;

typedef struct quadtree
{
  /* center of mass */
  double mass_tot;
  double cm_x;
  double cm_y;
  particle_t *particle;
  quadrant_t *quadrant;
  /* tree nodes */
  struct quadtree *top_left;
  struct quadtree *top_right;
  struct quadtree *bottom_left;
  struct quadtree *bottom_right;
}  quadtree_t;

int getInput(const char* __restrict filename, quadtree_t** __restrict root, particle_t* __restrict particles, const int N);
int insert(quadtree_t** __restrict node, particle_t* __restrict particle, quadrant_t* __restrict quadrant);
void insertQuadNode(quadtree_t** node, particle_t* __restrict particle);
void deleteQuadtree(quadtree_t **node);
void resetQuadtree(quadtree_t **node);
int setOutput(const char* __restrict filename, particle_t* __restrict data, const int N);
int calculateNewPositions(quadtree_t** __restrict root, particle_t* __restrict particles, const int N, const int nsteps, const double dt, const double theta_max);
void calculateForce(quadtree_t** __restrict node, particle_t* __restrict particle, double* __restrict F, const double theta_max);

/*
static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
*/

int main(int argc, char* argv[])
{
  //double time1 = get_wall_seconds();
  if (argc != 7)
  {
    printf("Check your input..\n");
    return -1;
  }
  const int N = atoi(argv[1]);
  const char* filename = argv[2];
  const int nsteps = atoi(argv[3]);
  const double dt = atof(argv[4]);
  const double theta_max = atof(argv[5]); // theta_max = 0.25
  //const int graphics = atoi(argv[6]);
  const char* filename_out = "result.gal";
  particle_t* particles = (particle_t*)malloc(N*sizeof(particle_t));
  quadtree_t* root = NULL;

  int err = getInput(filename, &root, particles, N);
  if(err == -1)
  {
    deleteQuadtree(&root); 
    free(particles);
    return -1;
  }
  err = calculateNewPositions(&root, particles, N, nsteps, dt, theta_max);
  if(err == -1)
  {
    deleteQuadtree(&root);
    free(particles);
    return -1;
  }
  
  err = setOutput(filename_out, particles, N);
  if(err == -1)
  {
    deleteQuadtree(&root);
    free(particles);
    return -1;
  }
  
  deleteQuadtree(&root);
  free(particles);
  //printf("galsim main took %7.3f wall seconds.\n", get_wall_seconds()-time1);
  return 0;
}

int getInput(const char* __restrict filename, quadtree_t** __restrict root, particle_t* __restrict particles, const int N)
{
  FILE *stream_in;
  stream_in = fopen(filename,"rb");
  double* arr = (double*)malloc(N*6*sizeof(double));
  quadrant_t *quadrant = (quadrant_t*)malloc(sizeof(particle_t));

  if(stream_in == NULL)
  {
    printf("Error: unable to open file: %s\n", filename);
    fclose(stream_in);
    free(arr);
    free(quadrant);
    return -1;
  }
  size_t input_size = N*6*sizeof(double);

  /* Read input to data array */
  size_t items_read = fread(arr, sizeof(char), input_size, stream_in);
  int j = 0;
  for (int i = 0; i < N*6; i += 6)
  {
    particles[j].pos_x = arr[i];
    particles[j].pos_x = arr[i];
    particles[j].pos_y = arr[i+1];
    particles[j].mass = arr[i+2];
    particles[j].vel_x = arr[i+3];
    particles[j].vel_y = arr[i+4];
    particles[j].brightness = arr[i+5];
    
    quadrant->width = 1;
    quadrant->height = 1;
    quadrant->origin_x = 0;
    quadrant->origin_y = 0;

    insert(root, &particles[j], quadrant);
    j++;
  }

  if (items_read != input_size)
  {
    printf("Error reading the input file.\n");
    fclose(stream_in);
    free(arr);
    free(quadrant);
    return -1;
  }
  free(quadrant);
  fclose(stream_in);
  free(arr);
  return 1;
}

/* Note: we define origin in bottom-left corner! */
int insert(quadtree_t** node, particle_t* __restrict particle, quadrant_t* __restrict quadrant)
{
  if(*node == NULL)
  {
    quadtree_t * new_node = (quadtree_t*)malloc(sizeof(quadtree_t));
    quadrant_t * new_quadrant = (quadrant_t*)malloc(sizeof(quadrant_t));
    new_node->particle = particle;
    new_node->mass_tot = particle->mass;
    new_node->cm_x = particle->pos_x;
    new_node->cm_y = particle->pos_y;

    new_quadrant->width = quadrant->width;
    new_quadrant->height = quadrant->height;
    new_quadrant->origin_x = quadrant->origin_x;
    new_quadrant->origin_y = quadrant->origin_y;
    new_node->quadrant = new_quadrant;

    new_node->top_left = NULL;
    new_node->top_right = NULL;
    new_node->bottom_left = NULL;
    new_node->bottom_right = NULL;
    *node = new_node;
  }
  /* 
   * internal node: a node which already contains other nodes 
   * external node: a node which is initially empty - after its been filled it will only contian one particle
   */
  else if ((*node)->particle == NULL) // internal node
  {
    /* updating center of mass */
    (*node)->cm_x = ((*node)->cm_x*(*node)->mass_tot + particle->pos_x*particle->mass)/((*node)->mass_tot + particle->mass);
    (*node)->cm_y = ((*node)->cm_y*(*node)->mass_tot + particle->pos_y*particle->mass)/((*node)->mass_tot + particle->mass);
    (*node)->mass_tot += particle->mass;
    insertQuadNode(&(*node), particle);
  }
  else // external node: move the new and old nodes into a new lower quadrant
  {
    if((*node)->particle->pos_x == particle->pos_x && (*node)->particle->pos_y == particle->pos_y){
      printf("Error: Particles at the same position \n");
      return -1;
    }
    /* Making our node internal by removing the particle */
    particle_t *particle_old = (*node)->particle;
    (*node)->particle = NULL;
    insertQuadNode(&(*node), particle);
    insertQuadNode(&(*node), particle_old);

    /* updating center of mass */
    (*node)->cm_x = ((*node)->cm_x*(*node)->mass_tot + particle->pos_x*particle->mass)/((*node)->mass_tot + particle->mass);
    (*node)->cm_y = ((*node)->cm_y*(*node)->mass_tot + particle->pos_y*particle->mass)/((*node)->mass_tot + particle->mass);
    (*node)->mass_tot += particle->mass;
  }
  return 0;
}

void insertQuadNode(quadtree_t** __restrict node, particle_t* __restrict particle)
{
  quadrant_t *quadrant_new = (quadrant_t*)malloc(sizeof(quadrant_t));
  quadrant_new->width = (*node)->quadrant->width/2.0;
  quadrant_new->height = (*node)->quadrant->height/2.0;

  if(particle->pos_x < ((*node)->quadrant->origin_x + quadrant_new->width)) // left quadrant
    {
      if (particle->pos_y < ((*node)->quadrant->origin_y + quadrant_new->height)) // bottom left quadrant
      {
        quadrant_new->origin_y = (*node)->quadrant->origin_y;
        quadrant_new->origin_x = (*node)->quadrant->origin_x;
        insert(&(*node)->bottom_left, particle, quadrant_new);
      }
      else // top left quadrant
      {
        quadrant_new->origin_y = (*node)->quadrant->origin_y + quadrant_new->height;
        quadrant_new->origin_x = (*node)->quadrant->origin_x;
        insert(&(*node)->top_left, particle, quadrant_new);
      }
    }
    else {  // right quadrant
      if (particle->pos_y < ((*node)->quadrant->origin_y + quadrant_new->height)) // bottom right quadrant
      {
        quadrant_new->origin_x = (*node)->quadrant->origin_x + quadrant_new->width;
        quadrant_new->origin_y = (*node)->quadrant->origin_y;
        insert(&(*node)->bottom_right, particle, quadrant_new);
      }
      else // top right quadrant
      {
        quadrant_new->origin_x = (*node)->quadrant->origin_x + quadrant_new->width;
        quadrant_new->origin_y = (*node)->quadrant->origin_y + quadrant_new->height;
        insert(&(*node)->top_right, particle, quadrant_new);
      }
    }
    free(quadrant_new);
}

void deleteQuadtree(quadtree_t **node)
{
   if(*node == NULL){ // if empty node
      return;
   }
   
  deleteQuadtree(&(*node)->top_left);
  deleteQuadtree(&(*node)->top_right);
  deleteQuadtree(&(*node)->bottom_left);
  deleteQuadtree(&(*node)->bottom_right);

  free((*node)->quadrant);
  free(*node);
  *node = NULL;
}

int setOutput(const char* __restrict filename, particle_t* __restrict data, const int N)
{
  FILE *stream_out;
  stream_out = fopen(filename,"wb");
  if(stream_out == NULL)
  {
    printf("Error: unable to open file: %s\n", filename);
    fclose(stream_out);
    return -1;
  }
  fwrite(data, sizeof(particle_t), N, stream_out);
  fclose(stream_out);
  return 0;
}

int calculateNewPositions(quadtree_t** root, particle_t* __restrict particles, const int N, const int nsteps, const double dt, const double theta_max)
{
  int err = 0;
  double ax, ay;
  const double G = 100.0/(double)(N);
  double* F = (double*)malloc(2*sizeof(double));
  particle_t *particles_new = (particle_t*)malloc(N*sizeof(particle_t));
  quadrant_t *quadrant = (quadrant_t*)malloc(sizeof(quadrant_t));
  quadtree_t **root_new = (quadtree_t**)malloc(sizeof(quadtree_t));
  *root_new = *root;

  for (int n = 0; n < nsteps; n++)
  {
    for(int i = 0; i < N; i++)
    {
      F[0] = 0;
      F[1] = 0;

      calculateForce(root_new, &particles[i], F, theta_max);

      F[0] = -G*particles[i].mass*F[0];
      F[1] = -G*particles[i].mass*F[1];
      
      /* Acceleration */
      ax = F[0]/(particles[i].mass);
      ay = F[1]/(particles[i].mass);
      /* Updating velocities */
      particles_new[i].vel_x = particles[i].vel_x + dt*ax;
      particles_new[i].vel_y = particles[i].vel_y + dt*ay;
      /* Updating positions */
      particles_new[i].pos_x = particles[i].pos_x + dt* particles_new[i].vel_x;
      particles_new[i].pos_y = particles[i].pos_y + dt* particles_new[i].vel_y;
      if(particles_new[i].pos_x >1|| particles_new[i].pos_x <0 || particles_new[i].pos_y >1|| particles_new[i].pos_y <0){
        printf("Error: Particle outside interval \n");
        free(particles_new);
        free(F);
        free(root_new);
        free(quadrant);
        return -1;
      }
    }
    deleteQuadtree(root_new);

    for (int i = 0; i < N; i++)
    {
      particles[i].pos_x = particles_new[i].pos_x;
      particles[i].pos_y = particles_new[i].pos_y;
      particles[i].vel_x = particles_new[i].vel_x;
      particles[i].vel_y = particles_new[i].vel_y;
      
      quadrant->width = 1;
      quadrant->height = 1;
      quadrant->origin_x = 0;
      quadrant->origin_y = 0;
      err += insert(root_new, &particles[i], quadrant);
    }
    if(err!=0)
    {
      free(particles_new);
      free(F);
      free(root_new);
      free(quadrant);
      return -1;
    }
  }
  *root = *root_new;
  free(particles_new);
  free(root_new);
  free(F);
  free(quadrant);
  return 0;
}

void calculateForce(quadtree_t** __restrict node, particle_t* __restrict particle, double* __restrict F, const double theta_max)
{
  double r, theta;
  const double epsilon = 0.001;
  /* F[0] = Fx, F[1] = Fy */
  if ((*node) == NULL)
  {
    return;
  }
  /* External node */
  else if ((*node)->particle != NULL)
  {
    /* Same particle */
    if ((fabs((*node)->particle->pos_x - particle->pos_x) < 0.0000001) && (fabs((*node)->particle->pos_y - particle->pos_y) < 0.0000001))
    {
      return;
    }
    /* x_i: particle, x_j: node->particle */
    r = sqrt((particle->pos_x - (*node)->particle->pos_x)*(particle->pos_x - (*node)->particle->pos_x) + (particle->pos_y - (*node)->particle->pos_y)*(particle->pos_y - (*node)->particle->pos_y));
    F[0] += (*node)->particle->mass*(particle->pos_x - (*node)->particle->pos_x)/((r+epsilon)*(r+epsilon)*(r+epsilon));
    F[1] += (*node)->particle->mass*(particle->pos_y - (*node)->particle->pos_y)/((r+epsilon)*(r+epsilon)*(r+epsilon));
  }
  /* Internal node */
  else
  {
    /* Calculate ratio */
    theta = (*node)->quadrant->width/sqrt((particle->pos_x - ((*node)->quadrant->origin_x + (*node)->quadrant->width/2))*(particle->pos_x - ((*node)->quadrant->origin_x + (*node)->quadrant->width/2)) + (particle->pos_y - ((*node)->quadrant->origin_y + (*node)->quadrant->height/2))*(particle->pos_y - ((*node)->quadrant->origin_y + (*node)->quadrant->height/2))); 
    
    /* Single body problem */
    if (theta <= theta_max) 
    {
      /* x_i: particle, x_j: node->cm */
      r = sqrt((particle->pos_x - (*node)->cm_x)*(particle->pos_x - (*node)->cm_x) + (particle->pos_y - (*node)->cm_y)*(particle->pos_y - (*node)->cm_y));
      F[0] += (*node)->mass_tot*(particle->pos_x - (*node)->cm_x)/((r+epsilon)*(r+epsilon)*(r+epsilon));
      F[1] += (*node)->mass_tot*(particle->pos_y - (*node)->cm_y)/((r+epsilon)*(r+epsilon)*(r+epsilon));
    }
    /* Many body problem */
    else
    {
      calculateForce(&(*node)->top_left, particle, F, theta_max);
      calculateForce(&(*node)->top_right, particle, F, theta_max);
      calculateForce(&(*node)->bottom_left, particle, F, theta_max);
      calculateForce(&(*node)->bottom_right, particle, F, theta_max);
    }
  }
} 
