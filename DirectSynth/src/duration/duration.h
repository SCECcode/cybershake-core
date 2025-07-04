#ifndef DURATION_H
#define DURATION_H

#define ARIAS_INTENSITY 0
#define ENERGY_INTEGRAL 1
#define CAV 2
#define DV 3
#define DA 4
#define D5_75 5
#define D5_95 6
#define D20_80 7

#define X_COMP 0
#define Y_COMP 1
#define Z_COMP 2

#define NUM_DURATION_MEASURES 9

#define NUM_DURATION_PERIODS 27

#define NUM_PERIOD_DURATION_MEASURES 3

struct duration_record {
  int type;
  int type_value;
  int component;
  float value;
};

struct period_duration_record {
	int type;
	int type_value;
	int component;
	float period;
	float value;
};

#endif
