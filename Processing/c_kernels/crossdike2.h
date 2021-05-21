struct point {
  double x;
  double y;
};

struct dike{
  struct point start;
  struct point stop;
};


bool inarray(int val, int arraysize, int *array){
  int i;
  for (i=0; i<arraysize; i++){
    if (val==*array) return true;
    array++;
  }
  return false;
}

int orientation(struct point a, struct point b, struct point c){
  double val = ((b.y - a.y) * (c.x - b.x)) - ((b.x - a.x) * (c.y - b.y));
  if (val == 0) return 0;  // colinear 
  return (val > 0)? 1: -1; // clock or counterclock wise
}

bool onSegment(struct point p, struct point q, struct point r) { 
    if (q.x <= fmax(p.x, r.x) && q.x >= fmin(p.x, r.x) && q.y <= fmax(p.y, r.y) && q.y >= fmin(p.y, r.y)) return true; 
    return false; 
}

bool intersect(struct point a1, struct point a2, struct point b1, struct point b2){
  // Function that returns true if line segment 'a' and 'b' intersect
  // Find the four orientations needed for general and special cases 
  int o1 = orientation(a1, a2, b1); 
  int o2 = orientation(a1, a2, b2); 
  int o3 = orientation(b1, b2, a1); 
  int o4 = orientation(b1, b2, a2); 
  // General case 
  if (o1 != o2 && o3 != o4)             return true; 
  // Special Cases 
  // a1, a2 and b1 are colinear and b1 lies on segment a1a2 
  if (o1 == 0 && onSegment(a1, b1, a2)) return true; 
  // a1, a2 and b2 are colinear and b2 lies on segment a1a2 
  if (o2 == 0 && onSegment(a1, b2, a2)) return true; 
  // b1, b2 and a1 are colinear and a1 lies on segment b1b2 
  if (o3 == 0 && onSegment(b1, a1, b2)) return true; 
  // b1, b2 and a2 are colinear and a2 lies on segment b1b2 
  if (o4 == 0 && onSegment(b1, a2, b2)) return true; 
  // Doesn't correspond to any of the above cases
  return false;
}

bool inside(int n, struct point polygon[], struct point p)
{ 
    // There must be at least 3 vertices in polygon[] 
    if (n < 3)  return false; 
  
    // Create a point for line segment from p to infinite (> 360° in earth coordonates)
    struct point extreme = {361, p.y}; 
  
    // Count intersections of the above line with sides of polygon 
    int count = 0, i = 0; 
    do
    { 
        int next = (i+1)%n; 
  
        // Check if the line segment from 'p' to 'extreme' intersects with the line segment from 'polygon[i]' to 'polygon[next]' 
        if (intersect(polygon[i], polygon[next], p, extreme)) 
        { 
            // If the point 'p' is colinear with line segment 'i-next', then check if it lies on segment. If it lies, return true, otherwise false 
            if (orientation(polygon[i], p, polygon[next]) == 0) 
               return onSegment(polygon[i], p, polygon[next]); 
            count++; 
        } 
        i = next; 
    } while (i != 0); 
  
    // Return true if count is odd, false otherwise 
    return count&1;  // Same as (count%2 == 1) 
} 

static inline StatusCode func(double *lon, double *lat, double *dep, int *beached,   double *prevlon, double *prevlat, double *prevdep)
                                                                               //int *points, double *dike_lon, double *dike_lat, int *restarts, int *restart)
{
  int p;
  struct point particle_start,particle_stop,dike_start,dike_stop;

  // in i,j coordinates (as seen in EFORIE tmask)
  // i=27 28 28 29 31 31 39 39 40 40 41 41   34  35  35  36     16  17  18    24 26
  // j=82 81 77 77 75 74 66 64 63 62 61 60  178 177 176 175    119 120 119    99 97
  int points=21;
  double dike_lon[21]={28.671237945556641,28.673706054687500,28.673706054687500,28.676176071166992,28.681114196777344,28.681114196777344,28.700866699218750,28.700866699218750,28.703336715698242,28.703336715698242,28.705804824829102,28.705804824829102,28.688522338867188,28.690990447998047,28.690990447998047,28.693460464477539,28.644077301025391,28.646545410156250,28.649015426635742,28.663829803466797,28.668767929077148};
  double dike_lat[21]={44.150741577148438,44.148887634277344,44.141479492187500,44.141479492187500,44.137775421142578,44.135925292968750,44.121109008789062,44.117404937744141,44.115554809570312,44.113704681396484,44.111850738525391,44.110000610351562,44.328517913818359,44.326667785644531,44.324813842773438,44.322963714599609,44.219257354736328,44.221111297607422,44.219257354736328,44.182220458984375,44.178516387939453};
  int restarts=3;
  int restart[3]={12,16,19};
  
  particle_start.x=*prevlon; particle_stop.x=*lon;
  particle_start.y=*prevlat; particle_stop.y=*lat;
  for (p=0; p<points; p++)  {
    if (restarts==0 || !(inarray(p+1,restarts,restart))){  // restart is equivalent to &restart[0]
      dike_start.x=dike_lon[p]; dike_stop.x=dike_lon[p+1];
      dike_start.y=dike_lat[p]; dike_stop.y=dike_lat[p+1];
      if (abs(*prevlon)>0.0000001)      {
	if (intersect(particle_start,particle_stop,dike_start,dike_stop)) {
	  //printf("dike collision %f %f -> %f %f\n",*prevlon,*prevlat,*lon,*lat);
          *lon = *prevlon;
	  *lat = *prevlat;
	  *dep = *prevdep;
	  break;
        }
      }
    }
  }

  // dike crossing check is last kernel component, so
  // save position as "previous position" to be used during next time step
  //*prevlon=*lon ;
  //*prevlat=*lat ;    !!! this is now a separate (python) kernel called after all the other ones
  //*prevdep=*dep ;

  return SUCCESS;
}
