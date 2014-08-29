#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <stdlib.h>
#include <math.h>

// simple 3D vector struct
struct vec3 {
  double x, y, z;
};

// simple but possibly inefficient getline (different from GNU getline)
int _getline(int fd, char *line, int size) {
  char c;
  int read_size = 0;
  while(size-- && read(fd, &c, 1)) {
    if(c == '\n') {
      *line = 0;
      read_size++;
      break;
    }
    *line = c;
    read_size++;
    line++;
  }
  return read_size;
}

// reads a vec3 from a CSV line
struct vec3 parse_vec3(char *line) {
  char *lasts;
  struct vec3 v;
  v.x = atof(strtok_r(line, ",", &lasts));
  v.y = atof(strtok_r(NULL, ",", &lasts));
  v.z = atof(lasts);
  return v;
}

/* --- some vector math functions --- */

struct vec3 sum(struct vec3 *v, int n) {
  struct vec3 s = {};
  while(n--) {
    s.x += v->x;
    s.y += v->y;
    s.z += v->z;
    v++;
  }
  return s;
}

double mag(struct vec3 v) {
  return sqrt(v.x * v.x +
	      v.y * v.y +
	      v.z * v.z);
}

struct vec3 scale(struct vec3 v, double a) {
  struct vec3 y = {
    .x = v.x * a,
    .y = v.y * a,
    .z = v.z * a
  };
  return y;
}

double rms(double *data, int size) {
  double ss = 0;
  int n = size;
  while(n--) {
    double x = *data++;
    ss += x * x;
  }
  return sqrt(ss / size);
}

double dot(struct vec3 a, struct vec3 b) {
  return
    a.x * b.x +
    a.y * b.y +
    a.z * b.z;
}

/* --- filtering functions --- */

// implements transpose direct form II filtering
double tdf2(unsigned int ord, const double *a, const double *b, double *z, const double x) {
  if(ord == 0) return b[0] * x;

  const double y = b[0] * x + z[0];
  for(int i = 1; i < ord; i++) {
    z[i-1] = b[i] * x + z[i] - a[i] * y;
  }
  z[ord-1] = b[ord] * x - a[ord] * y;
  return y;
}

// applies tdf2 on an array of data using calculated coefficients to implement
// a 4th order Butterworth bandpass filter for 1-3 Hz, given a 20 Hz sample rate.
void filter(double *x, int size, double *y) {
  const double b[9] = {0.00482434, 0, -0.01929737, 0, 0.02894606, 0, -0.01929737, 0, 0.00482434};
  const double a[9] = {1, -5.41823139, 13.5293587, -20.31926512, 20.07119886, -13.34437166, 5.83210677, -1.53473005, 0.18737949};
  double z[8] = {-0.00482434, -0.00482434, 0.01447303, 0.01447303, -0.01447303, -0.01447303, 0.00482434, 0.00482434};
  for(int i = 0; i < size; i++) {
    *y++ = tdf2(8, a, b, z, *x++);
  } 
}

// a simple threshold based counting function for the filtered signal
int count_steps(double *x, int size, double hi, double lo) {
  enum {
    UP = 0,
    DOWN
  } state = UP;
  
  int n = size;
  int cnt = 0;

  while(n--) {
    double v = *x++;
    switch(state) {
    case UP:
      if(v > hi) state = DOWN;
      break;
    case DOWN:
      if(v < lo) {
	state = UP;
	cnt++;
      }
      break;
    }
  }

  return cnt;
}

#define MAX_LINE_SIZE 1024
#define MAX_VECTORS 2048

// the main function takes a single argument from the command line: the CSV file name
int main(int argc, char **argv) {
  char line[MAX_LINE_SIZE];
  struct vec3 data[MAX_VECTORS];
  double vert[MAX_VECTORS];
  double filtered[MAX_VECTORS];
  int data_len = 0;

  if(argc < 2) {
    puts("No file provided.");
    return -1;
  }

  int fd = open(argv[1], O_RDONLY);

  // throw away the header line
  _getline(fd, line, sizeof(line));

  // parse and load the vector data
  struct vec3 *ptr = data;
  while(_getline(fd, line, sizeof(line)) && data_len < MAX_VECTORS) {
    struct vec3 v = parse_vec3(line);
    //printf("x = %f, y = %f, z = %f\n", v.x, v.y, v.z);
    *ptr++ = v;
    data_len++;
  }
  
  printf("vectors read: %d\n", data_len);

  // calculate the gravity vector
  struct vec3 s = sum(data, data_len);
  double m = mag(s);
  //printf("mag = %f\n", m);
  struct vec3 g = scale(s, 1/m);
  printf("normalized gravity vector: %f %f %f\n", g.x, g.y, g.z);

  // reduce 3D data to 1D vertical acceleration
  for(int i = 0; i < data_len; i++) {
    vert[i] = dot(data[i], g);
    //printf("v = %f\n", vert[i]);
  }

  // bandpass filter 1-3 Hz
  filter(vert, data_len, filtered);

  //for(int i = 0; i < data_len; i++) printf("%f\n", filtered[i]);
  
  // calculate the thresholds
  double rms_val = rms(filtered, data_len);
  printf("rms: %f\n", rms_val);
  double threshold = rms_val * 0.5;

  // count the steps in the cleaned up signal
  int cnt = count_steps(filtered, data_len, threshold, -threshold);  
  printf("cnt: %d\n", cnt);

  return 0;
}
