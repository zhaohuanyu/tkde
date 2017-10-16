#include "pointline.h"

double convexRedundancy(deque<Point> &upper,deque<Point> &lower);
void convexInitialization(deque<Point> &upper, deque<Point> &lower, Point point1, Point point2, Point tpoint1, Point tpoint2);
int feasibleRegionUpdate(deque<Point> &convex_upper, deque<Point> &convex_lower, Point point1, Point point2);
int feasibleRegionUpdateAddForward(deque<Point> &convex_upper, deque<Point> &convex_lower, Point point1, Point point2);
void semiOpt(string fileName, string resultName, string connectedIndex,string connectedValue, double sigma);
int disConnAlgOneVect(Line &intialLineU,Line &intialLineL,deque<Point> &convex_upper, deque<Point> &convex_lower,deque<Point> &convex_changeL, deque<Point> &convex_changeU, long startIndex, deque<Point> initialStream,double sigma,int &upOrDown, long &timeEnd, int upOrDown1, double &lt);
void optimalProcessingFS(string fileName, string resultName, string connectionIndex, string connectionValue, double sigma);
