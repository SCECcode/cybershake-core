#define PARL_HEAD 1
#define PARL_BODY 2
#define PARL_TAIL 3

#define SEND_VEL_LEFT	2
#define SEND_VEL_RIGHT	0
#define SEND_TAU_LEFT	6
#define SEND_TAU_RIGHT	4

#define RECV_VEL_LEFT	SEND_VEL_RIGHT
#define RECV_VEL_RIGHT	SEND_VEL_LEFT
#define RECV_TAU_LEFT	SEND_TAU_RIGHT
#define RECV_TAU_RIGHT	SEND_TAU_LEFT

#define MAX_PROC	99999

/*
   nodeType code is integer string of 3 places of the form:

      <X><Z><Y>

   where each of <X>,<Z>,<Y> are 0 (minimal boundary), 1 (interior) or 2 (maximal boundary)

   For example, if the node contains the x=0,z=0,y=0 boundary, then
      
      nodeType = 000

   a node along the x=0,z=0 boundary, but y=interior would be
      
      nodeType = 001

   a node along the x=xmax,z=0 boundary, but y=interior would be
      
      nodeType = 201

   a node x=z=y=interior would be
      
      nodeType = 111

   There are 27 possible combinations, as defined below.
*/

#define NODE_MIN_BOUNDARY	0
#define NODE_INTERIOR		1
#define NODE_MAX_BOUNDARY	2

#define NODE_XMIN_ZMIN_YMIN	000
#define NODE_XINT_ZMIN_YMIN	100
#define NODE_XMAX_ZMIN_YMIN	200

#define NODE_XMIN_ZINT_YMIN	010
#define NODE_XINT_ZINT_YMIN	110
#define NODE_XMAX_ZINT_YMIN	210

#define NODE_XMIN_ZMAX_YMIN	020
#define NODE_XINT_ZMAX_YMIN	120
#define NODE_XMAX_ZMAX_YMIN	220

#define NODE_XMIN_ZMIN_YINT	001
#define NODE_XINT_ZMIN_YINT	101
#define NODE_XMAX_ZMIN_YINT	201

#define NODE_XMIN_ZINT_YINT	011
#define NODE_XINT_ZINT_YINT	111
#define NODE_XMAX_ZINT_YINT	211

#define NODE_XMIN_ZMAX_YINT	021
#define NODE_XINT_ZMAX_YINT	121
#define NODE_XMAX_ZMAX_YINT	221

#define NODE_XMIN_ZMIN_YMAX	002
#define NODE_XINT_ZMIN_YMAX	102
#define NODE_XMAX_ZMIN_YMAX	202

#define NODE_XMIN_ZINT_YMAX	012
#define NODE_XINT_ZINT_YMAX	112
#define NODE_XMAX_ZINT_YMAX	212

#define NODE_XMIN_ZMAX_YMAX	022
#define NODE_XINT_ZMAX_YMAX	122
#define NODE_XMAX_ZMAX_YMAX	222

