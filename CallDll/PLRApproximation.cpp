#include "framework.h"
#include "fs.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

FS::FS() {
	lowerBound = 0;
	upperBound = 0;
}

FS::FS(const FS & inObj) {
	lowerBound = inObj.lowerBound;
	upperBound = inObj.upperBound;
}

Point pointCreate(double xVal, double yVal)
{
	Point newPoint;
	newPoint.x = xVal;
	newPoint.y = yVal;
	return newPoint;
}

void pointCreate(Point &newPoint, double xVal, double yVal)
{
	newPoint.x=xVal;
	newPoint.y=yVal;
	return;
}

Line lineCreate(double xVal1, double yVal1, double xVal2, double yVal2)
{
	Line newLine;
	if(xVal1==xVal2)
	{
		newLine.slope = 0;
		newLine.offset = (yVal1+yVal2)/2;
	}
	else
	{
		newLine.slope = (yVal1 - yVal2)/(xVal1 - xVal2);
		newLine.offset= (xVal1*yVal2 - yVal1*xVal2)/(xVal1 - xVal2);
	}
	return newLine;
}

Line lineCreate(const Point &point1, const Point &point2)
{
	Line newLine;
	if(point1.x==point2.x)
	{
		newLine.slope=0;
		newLine.offset=(point1.y+point2.y)/2;
	}
	else
	{
		newLine.slope = (point1.y - point2.y)/(point1.x - point2.x);
		newLine.offset= (point1.x*point2.y - point1.y*point2.x)/(point1.x - point2.x);
	}
	return newLine;
}

Line lineCreate(double slopeVal, const Point &point)
{
	Line newLine;
	newLine.slope = slopeVal;
	newLine.offset = point.y - slopeVal * point.x;
	return newLine;
}

Line lineCreate(double slopeVal, double offsetVal)
{
	Line newLine;
	newLine.slope = slopeVal;
	newLine.offset = offsetVal;
	return newLine;
}

void lineCreate(Line &newLine, double xVal1, double yVal1, double xVal2, double yVal2)
{
	if(xVal1==xVal2)
	{
		newLine.slope = 0;
		newLine.offset = (yVal1+yVal2)/2;
	}
	else
	{
		newLine.slope = (yVal1 - yVal2)/(xVal1 - xVal2);
		newLine.offset= (xVal1*yVal2 - yVal1*xVal2)/(xVal1 - xVal2);
	}
	return;
}

void lineCreate(Line &newLine, const Point &point1, const Point &point2)
{
	if(point1.x==point2.x)
	{
		newLine.slope=0;
		newLine.offset=(point1.y+point2.y)/2;
	}
	else
	{
		newLine.slope = (point1.y - point2.y)/(point1.x - point2.x);
		newLine.offset= (point1.x*point2.y - point1.y*point2.x)/(point1.x - point2.x);
	}
	return;
}

void lineCreate(Line &newLine, double slopeVal, double offsetVal)
{
	newLine.slope = slopeVal;
	newLine.offset = offsetVal;
	return;
}

Point lineIntersection(const Line &line1, const Line &line2)
{
	Point interP;
	if(line1.slope==line2.slope)
	{
		interP.x=0;
		interP.y=line1.offset;
	}
	else
	{
		interP.x = (line1.offset - line2.offset)/(line2.slope - line1.slope);
		interP.y = (line2.slope*line1.offset - line1.slope*line2.offset)/(line2.slope - line1.slope);
	}
	return interP;
}

void lineIntersection(Point &interP, const Line &line1, const Line &line2)
{
	if(line1.slope==line2.slope)
	{
		interP.x=0;
		interP.y=line1.offset;
	}
	else
	{
		interP.x = (line1.offset - line2.offset)/(line2.slope - line1.slope);
		interP.y = (line2.slope*line1.offset - line1.slope*line2.offset)/(line2.slope - line1.slope);
	}
	return;
}

double linePoint(const Line &line, double xVal)
{
	double yVal;
	yVal = line.slope*xVal+line.offset;
	return yVal;
}

double getSlope(const Point &point1, const Point &point2)
{
	double slope;
	if(point1.x == point2.x)
	{
		slope = 0;
	}
	else
	{
		slope=(point1.y - point2.y)/(point1.x - point2.x);
	}
	return slope;
}

double convexRedundancy(deque<Point> &upper,deque<Point> &lower)
{
	//need revise
	double redundancyRate;
	double upperRedundantPointNumber=0;
	double lowerRedundantPointNumber=0;
	int i;
	
	if(upper.size()>=3)
	{
		for(i=2;i<upper.size();i++)
		{
			if(abs(getSlope(upper[i-2],upper[i-1])-getSlope(upper[i-1],upper[i])) < 0.000001)
			{
				upperRedundantPointNumber++;
			}
		}
	}
	if(lower.size()>=3)
	{
		for(i=2;i<lower.size();i++)
		{
			if(abs(getSlope(lower[i-2],lower[i-1])-getSlope(lower[i-1],lower[i])) < 0.000001)
			{
				lowerRedundantPointNumber++;
			}
		}
	}
	
	redundancyRate = (upperRedundantPointNumber+lowerRedundantPointNumber)/(upper.size()+lower.size());
	return redundancyRate;
}

void convexInitialization(deque<Point> &upper, deque<Point> &lower, Point point1, Point point2, Point tpoint1, Point tpoint2)
{
	upper.push_back(point1);upper.push_back(tpoint1);
	lower.push_back(point2);lower.push_back(tpoint2);
	return;
}

int feasibleRegionUpdate(deque<Point> &convex_upper, deque<Point> &convex_lower, Point point1, Point point2)
{
	//Line tempL, tempU;
	//tempL = lineCreate(convex_upper.front(), convex_lower.back());
	//tempU = lineCreate(convex_upper.back(), convex_lower.front());

	//in order to remove the effect of division operation, we use the ulternatives in if and else

	//update point 1
	//if(linePoint(tempL,point1.x) > point1.y)
	if((point1.x-convex_lower.back().x)*(convex_lower.back().y-convex_upper.front().y)>(point1.y-convex_lower.back().y)*(convex_lower.back().x-convex_upper.front().x))
	{
		return 0;
	}
	//else if(linePoint(tempU,point1.x) <= point1.y)
	else if((point1.x-convex_upper.back().x)*(convex_upper.back().y-convex_lower.front().y)<=(point1.y-convex_upper.back().y)*(convex_upper.back().x-convex_lower.front().x))
	{
		//do nothing
	}
	else
	{
		while(convex_lower.size()>1 && (convex_lower[1].y-convex_lower[0].y)*(point1.x-convex_lower[0].x)>=(point1.y-convex_lower[0].y)*(convex_lower[1].x-convex_lower[0].x))
		{
			convex_lower.pop_front();
		}
		
		while(convex_upper.size()>1 && (convex_upper[convex_upper.size()-1].y-convex_upper[convex_upper.size()-2].y)*(point1.x-convex_upper[convex_upper.size()-1].x)>=(point1.y-convex_upper[convex_upper.size()-1].y)*(convex_upper[convex_upper.size()-1].x-convex_upper[convex_upper.size()-2].x))
		{
			convex_upper.pop_back();
		}
		
		convex_upper.push_back(point1);
	}
	
	//update point 2
	//if(linePoint(tempU,point2.x) < point2.y)
	if((point2.x-convex_upper.back().x)*(convex_upper.back().y-convex_lower.front().y)<(point2.y-convex_upper.back().y)*(convex_upper.back().x-convex_lower.front().x))
	{
		return 0;
	}
	//else if(linePoint(tempL,point2.x) >= point2.y)
	else if((point2.x-convex_lower.back().x)*(convex_lower.back().y-convex_upper.front().y)>=(point2.y-convex_lower.back().y)*(convex_lower.back().x-convex_upper.front().x))
	{
		//do nothing
	}
	else
	{
		while(convex_upper.size()>1 && (convex_upper[1].y-convex_upper[0].y)*(point2.x-convex_upper[0].x)<=(point2.y-convex_upper[0].y)*(convex_upper[1].x-convex_upper[0].x))
		{
			convex_upper.pop_front();
		}
		
		while(convex_lower.size()>1 && (convex_lower[convex_lower.size()-1].y-convex_lower[convex_lower.size()-2].y)*(point2.x-convex_lower[convex_lower.size()-1].x)<=(point2.y-convex_lower[convex_lower.size()-1].y)*(convex_lower[convex_lower.size()-1].x-convex_lower[convex_lower.size()-2].x))
		{
			convex_lower.pop_back();
//			if(point2.x==67)
//			{
//				Point pt1 = convex_lower[convex_lower.size()-1];
//				Point pt2 = convex_lower[convex_lower.size()-2];
//				int aa =1;
//			}
		}

		convex_lower.push_back(point2);
	}
	
	return 1;
}

int feasibleRegionUpdateAddForward(deque<Point> &convex_upper, deque<Point> &convex_lower, Point point1, Point point2)
{
	Line tempL, tempU;
	//coutDque(convex_upper,convex_lower);
	tempL = lineCreate(convex_upper.front(), convex_lower.back());
	tempU = lineCreate(convex_upper.back(), convex_lower.front());

	//in order to remove the effect of division operation, we use the ulternatives in if and else

	//update point 1
	double dtemp1 = linePoint(tempU,point1.x);
	double dtemp2 = linePoint(tempL,point1.x);
	//if(linePoint(tempU,point1.x) > point1.y)
	if((convex_lower.front().x-point1.x)*(convex_upper.back().y-convex_lower.front().y)<(convex_lower.front().y-point1.y)*(convex_upper.back().x-convex_lower.front().x))
	{
		return 0;
	}
	//else if(linePoint(tempL,point1.x) <= point1.y)
	else if((convex_upper.front().x-point1.x)*(convex_lower.back().y-convex_upper.front().y)>=(convex_upper.front().y-point1.y)*(convex_lower.back().x-convex_upper.front().x))
	{
		//do nothing
	}
	else
	{
		//把convex hull下界最右边的点后面的点删除
		while(convex_lower.size()>1 && (convex_lower[convex_lower.size()-1].y-convex_lower[convex_lower.size()-2].y)*(convex_lower[convex_lower.size()-1].x-point1.x)<=(convex_lower[convex_lower.size()-1].y-point1.y)*(convex_lower[convex_lower.size()-1].x-convex_lower[convex_lower.size()-2].x))
		{
			convex_lower.pop_back();
		}
		//三角检验convex hull的上界左边的点
		while(convex_upper.size()>1 && (convex_upper[1].y-convex_upper[0].y)*(convex_upper[0].x-point1.x)<=(convex_upper[0].y-point1.y)*(convex_upper[1].x-convex_upper[0].x))
		{
			convex_upper.pop_front();
		}

		convex_upper.push_front(point1);
	}

	//update point 2
	dtemp1 = linePoint(tempL,point2.x);
	dtemp2 = linePoint(tempU,point2.x);
	//if(linePoint(tempL,point2.x) < point2.y)
	if((convex_upper.front().x-point2.x)*(convex_lower.back().y-convex_upper.front().y)>(convex_upper.front().y-point2.y)*(convex_lower.back().x-convex_upper.front().x))
	{
		return 0;
	}
	//else if(linePoint(tempU,point2.x) >= point2.y)
	else if((convex_lower.front().x-point2.x)*(convex_upper.back().y-convex_lower.front().y)<=(convex_lower.front().y-point2.y)*(convex_upper.back().x-convex_lower.front().x))
	{
		//do nothing
	}
	else
	{
		while(convex_upper.size()>1 && (convex_upper[convex_upper.size()-1].y-point2.y)*(convex_upper[convex_upper.size()-1].x-convex_upper[convex_upper.size()-2].x)<=(convex_upper[convex_upper.size()-1].x-point2.x)*(convex_upper[convex_upper.size()-1].y-convex_upper[convex_upper.size()-2].y))
		{
			convex_upper.pop_back();
		}

		while(convex_lower.size()>1 && (convex_lower[1].y-convex_lower[0].y)*(convex_lower[0].x-point2.x)>=(convex_lower[0].y-point2.y)*(convex_lower[1].x-convex_lower[0].x))
		{
			convex_lower.pop_front();
		}

		convex_lower.push_front(point2);
	}

	//coutDque(convex_upper,convex_lower);
	tempL = lineCreate(convex_upper.front(), convex_lower.back());
	tempU = lineCreate(convex_upper.back(), convex_lower.front());
	return 1;
}

int disConnAlgOneVect(Line &intialLineU,Line &intialLineL,deque<Point> &convex_upper, deque<Point> &convex_lower,deque<Point> &convex_changeL, deque<Point> &convex_changeU, long startIndex, deque<Point> initialStream,double sigma,int &upOrDown, long &timeEnd, int upOrDown1)
{
	double streamdata;
	int sIndexOffset=0;
	long lindex = startIndex;
	Point point1, point2, tpoint1, tpoint2;
	convex_changeL.clear();
	convex_changeU.clear();

	while(startIndex<initialStream.size())
	{
		sIndexOffset++;
		streamdata=initialStream[startIndex-1].y;
		pointCreate(point1,startIndex,streamdata+sigma);//upper point
		pointCreate(point2,startIndex,streamdata-sigma);//lower point
		if(sIndexOffset == 1)
		{
			if(startIndex<initialStream.size())
			{
				startIndex++;
				sIndexOffset++;
				streamdata=initialStream[startIndex-1].y;
				pointCreate(tpoint1,startIndex,streamdata+sigma);//upper point
				pointCreate(tpoint2,startIndex,streamdata-sigma);//lower point
				convexInitialization(convex_upper,convex_lower,point1,point2,tpoint1,tpoint2);
				lineCreate(intialLineL,convex_upper.back(), convex_lower.front());
				lineCreate(intialLineU,convex_upper.front(), convex_lower.back());
				//初始化记录convex_hull变化的位置
				Point ptU = pointCreate(convex_upper.front().x,convex_lower.back().x);
				Point ptL = pointCreate(convex_lower.front().x,convex_upper.back().x);
				convex_changeL.push_back(ptL);
				convex_changeU.push_back(ptU);
				startIndex++;
				continue;
			}
			else
			{
				cout<<"Note the last segment contains only one point!"<<endl;
				return 0;
			}
		}

		if(!feasibleRegionUpdate(convex_upper,convex_lower,point1,point2))
		{
			//输出方向
			Line tempL, tempU;//tempL:上界；tempU：下界
			lineCreate(tempL,convex_upper.back(), convex_lower.front());
			lineCreate(tempU,convex_upper.front(), convex_lower.back());
			//下界的交点在tolerate+上，是从下往上
			if(tempU.slope * point1.x + tempU.offset > point1.y)
			{
				upOrDown=0;
			}
			//上界的交点在tolerate-下，是从上往下
			if(tempL.slope * point1.x + tempL.offset < point2.y)
			{
				upOrDown=1;
			}
			timeEnd=point1.x;
			return 1;
		}
		else
		{
			//if(lindex==61)
			//coutDque(convex_upper,convex_lower);
			if(upOrDown1 == 1)
			{
				Point ptTemp1 = convex_changeL.back();
				Point ptTemp2 = convex_upper.back();
				if(ptTemp1.y < ptTemp2.x)
				{
					Point ptL = pointCreate(convex_lower.front().x,convex_upper.back().x);
					convex_changeL.push_back(ptL);
				}
			}
			if(upOrDown1 == 0)
			{
				Point ptTemp1 = convex_changeU.back();
				Point ptTemp2 = convex_lower.back();
				if(ptTemp1.y < ptTemp2.x)
				{
					Point ptU = pointCreate(convex_upper.front().x,convex_lower.back().x);
					convex_changeU.push_back(ptU);
				}
			}

			startIndex++;
			if(startIndex<initialStream.size())
			{
				continue;
			}
			else
			{
				cout<<"Note the last segment contains only one point!"<<endl;
				return 0;
			}
		}
	}
}

void semiOpt(string fileName, string resultName, string connectionIndex, string connectionValue, double sigma)
{

	ifstream streamFile(fileName.c_str());
	ofstream resultFile(resultName.c_str(),ios::app);
	ofstream connectIndexFile(connectionIndex.c_str());
	ofstream connectFile(connectionValue.c_str());

	//first read all the sequence into memory
	deque<Point> initialStream;
	Point initialPoint;
	initialPoint=pointCreate(0,0);
	long initialIndex=1;
	double initialData;
	double maxVal, minVal;

	maxVal = -1000000000;
	minVal = 10000000000;

	FILE *fp;
	fp=fopen(fileName.c_str(),"r");
	if (fp == NULL)
	{
	    exit(1);
	}
	while(!feof(fp))
	{
		try
		{
		 int f = fscanf(fp,"%lf",&initialData);
	     pointCreate(initialPoint,initialIndex,initialData);
	     if(f == 0)
	     {
	    	 exit(1);
	     }
	     if(maxVal < initialData)
	     {
	     	maxVal = initialData;
	     }
	     if(minVal > initialData)
	     {
	     	minVal = initialData;
	     }
	     initialStream.push_back(initialPoint);
	     initialIndex++;
		}
			catch(exception e)
			{
				exit(1);
			}
	}
	fclose(fp);
	resultFile<<"误差："<< " "<< sigma << endl;
	resultFile<<""<< " " << endl;
	deque<Point> convex_upper, convex_lower; //convex_upper is located on the above
	deque<Point> convex_changeL, convex_changeU;
	Point point1, point2, tpoint1, tpoint2;
	Line tempU, tempL; //record the segment slope with maximum value
	double streamdata;
	long sIndex; //record the index of the point under processing
	long startIndex; //record the index of the first point of a segment
	long sIndexOffset; //record the relative index of the point in current segment
	long numberSeg; //record the number of segments

	tempU=lineCreate(0,0);
	streamdata=0;
	sIndex = 0;
	startIndex = 1;
	sIndexOffset = 0;
	numberSeg = 0;

	clock_t start,finish, startTest, finishTest, startTest1, finishTest1;
	double duration;
	start = clock();
	sIndex=0;

	int upOrDown,currentUpOrDown;
	long timeEnd;
	int aa=0;
	deque<Point> connectPoint;
	Line keyLine, currentIntialLineU, currentIntialLineL,nextIntialLineU, nextIntialLineL;
	Point keyPoint;
	double lTimeTest=0;
	double ltimeTest1=0;

	//初始化第一个分割
	disConnAlgOneVect(currentIntialLineU,currentIntialLineL,convex_upper, convex_lower, convex_changeL,convex_changeU,startIndex, initialStream, sigma, upOrDown, timeEnd, 2);
	//记录key extremes line和key convex hull point
	if(upOrDown == 1)
	{
		//back指的是靠后，front指的是考前
		 lineCreate(keyLine,convex_lower.front(), convex_upper.back());
		 keyPoint = convex_upper.back();
		 connectPoint.push_back(pointCreate(startIndex,keyLine.slope * startIndex + keyLine.offset));
		 numberSeg++;
		 if(numberSeg % 257 == 0)
		 {
			 connectFile << endl;
			 connectIndexFile << endl;
		 }
		 connectIndexFile << startIndex << "    ";
		 connectFile << keyLine.slope * startIndex + keyLine.offset << "    ";
	}
	if(upOrDown == 0)
	{
		 lineCreate(keyLine,convex_upper.front(), convex_lower.back());
		 keyPoint = convex_lower.back();
		 connectPoint.push_back(pointCreate(startIndex,keyLine.slope * startIndex + keyLine.offset));
		 numberSeg++;
		 if(numberSeg % 257 == 0)
		 {
			 connectFile << endl;
			 connectIndexFile << endl;
		 }
		 connectIndexFile << startIndex << "    ";
		 connectFile << keyLine.slope * startIndex + keyLine.offset << "    ";
	}
	startIndex=timeEnd;
	convex_lower.clear();
	convex_upper.clear();
	bool blwhile=true;
	bool lastseg=false;
	bool blonlyone=false;
	while(blwhile)
	{
		double streamdata;
		int sIndexOffset=0;
		long lindex = startIndex;
		bool blwhile1=true;
		Point point1, point2, tpoint1, tpoint2;
		convex_changeL.clear();
		convex_changeU.clear();
		convex_upper.clear();
		convex_lower.clear();
		while(lindex<=initialStream.size())
		{
			sIndexOffset++;
			streamdata=initialStream[lindex-1].y;
			pointCreate(point1,lindex,streamdata+sigma);//upper point
			pointCreate(point2,lindex,streamdata-sigma);//lower point
			if(sIndexOffset == 1)
			{
				if(lindex<initialStream.size())
				{
					lindex++;
					sIndexOffset++;
					streamdata=initialStream[lindex-1].y;
					pointCreate(tpoint1,lindex,streamdata+sigma);//upper point
					pointCreate(tpoint2,lindex,streamdata-sigma);//lower point
					convexInitialization(convex_upper,convex_lower,point1,point2,tpoint1,tpoint2);
					lineCreate(nextIntialLineL,convex_upper.back(), convex_lower.front());
					lineCreate(nextIntialLineU,convex_upper.front(), convex_lower.back());
					//初始化记录convex_hull变化的位置
					Point ptU = pointCreate(convex_upper.front().x,convex_lower.back().x);
					Point ptL = pointCreate(convex_lower.front().x,convex_upper.back().x);
					convex_changeL.push_back(ptL);
					convex_changeU.push_back(ptU);
					lindex++;
					if(lindex>initialStream.size())
					{
						blwhile1=true;
						lastseg=true;
						break;
					}
					continue;
				}
				else
				{
					cout<<"Note the last segment contains only one point!"<<endl;
					lindex++;
					blonlyone=true;
					break;
				}
			}
			if(!feasibleRegionUpdate(convex_upper,convex_lower,point1,point2))
			{
				//输出方向
				Line tempL, tempU;//tempL:上界；tempU：下界
				lineCreate(tempL,convex_upper.back(), convex_lower.front());
				lineCreate(tempU,convex_upper.front(), convex_lower.back());
				//下界的交点在tolerate+上，是从下往上
				if(tempU.slope * point1.x + tempU.offset > point1.y)
				{
					currentUpOrDown=0;
				}
				//上界的交点在tolerate-下，是从上往下
				if(tempL.slope * point1.x + tempL.offset < point2.y)
				{
					currentUpOrDown=1;
				}
				timeEnd=point1.x;
				blwhile1=true;
				break;
			}
			else
			{
				if(upOrDown == 1)
				{
					Point ptTemp1 = convex_changeL.back();
					Point ptTemp2 = convex_upper.back();
					if(ptTemp1.y < ptTemp2.x)
					{
						Point ptL = pointCreate(convex_lower.front().x,ptTemp2.x);
						convex_changeL.push_back(ptL);
					}
				}
				if(upOrDown == 0)
				{
					Point ptTemp1 = convex_changeU.back();
					Point ptTemp2 = convex_lower.back();
					if(ptTemp1.y < ptTemp2.x)
					{
						Point ptU = pointCreate(convex_upper.front().x,ptTemp2.x);
						convex_changeU.push_back(ptU);
					}
				}
				lindex++;
				if(lindex<=initialStream.size())
				{
					continue;
				}
				else
				{
					blwhile1=true;
					lastseg=true;
					break;
				}
			}
		}
		if(blonlyone)
		{
			connectPoint.push_back(pointCreate(initialStream.size()-1,keyLine.slope * (initialStream.size()-1) + keyLine.offset));
			numberSeg++;
			 if((numberSeg) % 257 == 0)
			 {
				 connectFile << endl;
				 connectIndexFile << endl;
			 }
			 connectIndexFile << initialStream.size()-1 << "    ";
			 connectFile << keyLine.slope * (initialStream.size()-1) + keyLine.offset << "    ";
			connectPoint.push_back(pointCreate(initialStream.size(),initialStream[initialStream.size()-1].y));
			numberSeg++;
			if(numberSeg % 257 == 0)
			 {
				 connectFile << endl;
				 connectIndexFile << endl;
			 }
			connectIndexFile << initialStream.size() << "    ";
			connectFile << initialStream[initialStream.size()-1].y << "    ";
			break;
		}

		//判断初始化的两条边界直线中靠近上一个分割key line的是否有交点
		//判断方法：先判断两条直线的交点是否大于等于上一个分割的keypoint的位置，如果小于，直接说明没有交点，如果大于等于，判断之前的是否在sigma范围内
		int isIntersect=1;//0表示没有交点，1表示有焦点
		Point pointIntersect;
		if(upOrDown == 1)
		{
			pointIntersect = lineIntersection(nextIntialLineL, keyLine);
			if(pointIntersect.x >= keyPoint.x && pointIntersect.x < startIndex&& abs(pointIntersect.x-startIndex) > 0.0001)
			{
				for(int i=startIndex-1; i > (int)(pointIntersect.x); i--)
				{
					if(nextIntialLineL.slope * i+nextIntialLineL.offset > initialStream[i-1].y+sigma)
					{
						isIntersect = 0;
						break;
					}
				}
			}
			else
			{
				isIntersect=0;
			}
		}

		if(upOrDown == 0)
		{
			pointIntersect = lineIntersection(nextIntialLineU, keyLine);
			if(pointIntersect.x >= keyPoint.x && pointIntersect.x < startIndex&& abs(pointIntersect.x-startIndex) > 0.0001)
			{
				for(int i=startIndex-1; i > (int)(pointIntersect.x); i--)
				{
					if(nextIntialLineU.slope * i+nextIntialLineU.offset < initialStream[i-1].y-sigma)
					{
						isIntersect = 0;
						break;
					}
				}
			}
			else
			{
				isIntersect=0;
			}
		}

		if(isIntersect == 0)
		{
			//不相交的情况下，直接保存，上一个分割的key extremes line在该分隔末点的交点
			//此外还要保存当前分割的key extremes line在初始点的交点
			 connectPoint.push_back(pointCreate(startIndex-1,keyLine.slope * (startIndex-1) + keyLine.offset));
			 numberSeg++;
			 if((numberSeg) % 257 == 0)
			 {
				 connectFile << endl;
				 connectIndexFile << endl;
			 }
			 connectIndexFile << startIndex-1<< "    ";
			 connectFile << keyLine.slope * (startIndex-1) + keyLine.offset << "    ";
			if(currentUpOrDown == 1)
			{
				 lineCreate(keyLine,convex_upper.back(), convex_lower.front());
				 keyPoint = convex_upper.back();
			}
			if(currentUpOrDown == 0)
			{
				 lineCreate(keyLine,convex_lower.back(), convex_upper.front());
				 keyPoint = convex_lower.back();
			}
			connectPoint.push_back(pointCreate(startIndex,keyLine.slope * startIndex + keyLine.offset));
			numberSeg++;
			 if((numberSeg) % 257 == 0)
			 {
				 connectFile << endl;
				 connectIndexFile << endl;
			 }
			 connectIndexFile << startIndex << "    ";
			 connectFile << keyLine.slope * (startIndex) + keyLine.offset << "    ";
			 if(lastseg==true)
			 {
				if(lindex>initialStream.size())
					timeEnd=lindex;
				connectPoint.push_back(pointCreate(lindex-1,keyLine.slope * (lindex-1) + keyLine.offset));
				numberSeg++;
				 if((numberSeg) % 257 == 0)
				 {
					 connectFile << endl;
					 connectIndexFile << endl;
				 }
				 connectIndexFile << lindex-1 << "    ";
				 connectFile << keyLine.slope * (lindex-1) + keyLine.offset << "    ";
				blwhile=false;
			 }
			startIndex=timeEnd;
			upOrDown=currentUpOrDown;
			convex_lower.clear();
			convex_upper.clear();
			continue;
		}

		if(isIntersect == 1)
		{
			//判断是否有交点
			int isLongIntersect=1;//0表示没有交点，1表示有焦点
			Point pointLongIntersect;
			Line currentLineU=lineCreate(convex_lower.back(), convex_upper.front());
			Line currentLineL=lineCreate(convex_upper.back(), convex_lower.front());
			if(upOrDown == 1)
			{
				pointLongIntersect = lineIntersection(currentLineL, keyLine);
				if(pointLongIntersect.x >= keyPoint.x && pointLongIntersect.x < startIndex)
				{
					for(int i=startIndex-1; i >(int)(pointLongIntersect.x); i--)
					{
						if(currentLineL.slope * i+currentLineL.offset > initialStream[i-1].y+sigma)
						{
							isLongIntersect = 0;
							break;
						}
					}
				}
				else
				{
					isLongIntersect=0;
					if(lastseg==true)
					{
						lastseg=false;
					}
				}
			}

			if(upOrDown == 0)
			{
				pointLongIntersect = lineIntersection(currentLineU, keyLine);
				if(pointLongIntersect.x >= keyPoint.x && pointLongIntersect.x < startIndex)
				{
					for(int i=startIndex-1; i > (int)(pointLongIntersect.x); i--)
					{
						if(currentLineU.slope * i+currentLineU.offset < initialStream[i-1].y-sigma)
						{
							isLongIntersect = 0;
							break;
						}
					}
				}
				else
				{
					isLongIntersect=0;
					if(lastseg==true)
					{
						lastseg=false;
					}
				}
			}

			if(isLongIntersect==0)
			{
				startTest1 = clock();
				//没有交点的话，从当前分割前到后去点直到有交点
				//注意，更新当前的convex hulls和key extremes line以及key convex hull
				Point pointTempS, pointTempX, pointEndS;
				int isLongIntersectTemp=0;
				int endIndex = timeEnd-1;
				Point pointLongIntersectTemp;
				Point convexChangeLPoint=convex_changeL.back();
				Point convexChangeUPoint=convex_changeU.back();
				Line currentLineUMax = lineCreate(pointCreate(convexChangeUPoint.x,initialStream[convexChangeUPoint.x-1].y+sigma),pointCreate(convexChangeUPoint.y,initialStream[convexChangeUPoint.y-1].y-sigma));
				Line currentLineLMax = lineCreate(pointCreate(convexChangeLPoint.x,initialStream[convexChangeLPoint.x-1].y-sigma),pointCreate(convexChangeLPoint.y,initialStream[convexChangeLPoint.y-1].y+sigma));
				convex_changeL.pop_back();
				convex_changeU.pop_back();
				while(isLongIntersectTemp==0)
				{
					if(upOrDown == 1)
					{
						Point convexChangeLPointTemp=convex_changeL.back();
						Line currentLineLTemp=lineCreate(pointCreate(convexChangeLPointTemp.x,initialStream[convexChangeLPointTemp.x-1].y-sigma),pointCreate(convexChangeLPointTemp.y,initialStream[convexChangeLPointTemp.y-1].y+sigma));
							currentLineLMax.slope=currentLineLTemp.slope;
							currentLineLMax.offset=currentLineLTemp.offset;
							pointLongIntersectTemp = lineIntersection(currentLineLTemp, keyLine);
							if(pointLongIntersectTemp.x >= keyPoint.x && pointLongIntersectTemp.x < startIndex)
							{
								bool bTemp=true;
								for(int i=startIndex; i > (int)pointLongIntersectTemp.x; i--)
								{
									double dbtmpe = initialStream[i-1].y+sigma;
									if(currentLineLTemp.slope * i+currentLineLTemp.offset > initialStream[i-1].y+sigma)
									{
										bTemp=false;
										break;
									}
								}
								if(bTemp)
								{
									isLongIntersectTemp=1;
									endIndex = convexChangeLPointTemp.y;
								}
							}
					}

					if(upOrDown == 0)
					{
						Point convexChangeUPointTemp=convex_changeU.back();
						Line currentLineUTemp=lineCreate(pointCreate(convexChangeUPointTemp.x,initialStream[convexChangeUPointTemp.x-1].y+sigma),pointCreate(convexChangeUPointTemp.y,initialStream[convexChangeUPointTemp.y-1].y-sigma));
							currentLineUMax.slope=currentLineUTemp.slope;
							currentLineUMax.offset=currentLineUTemp.offset;
							pointLongIntersectTemp = lineIntersection(currentLineUTemp, keyLine);
							if(pointLongIntersectTemp.x >= keyPoint.x && pointLongIntersectTemp.x < startIndex)
							{
								bool bTemp=true;
								for(int i=startIndex; i > (int)pointLongIntersectTemp.x; i--)
								{
									double dbtmpe = initialStream[i-1].y+sigma;
									if(currentLineUTemp.slope * i+currentLineUTemp.offset < initialStream[i-1].y-sigma)
									{
										bTemp=false;
										break;
									}
								}
								if(bTemp)
								{
									isLongIntersectTemp=1;
									endIndex = convexChangeUPointTemp.y;
								}
							}
						//}
					}
					if(isLongIntersectTemp==0)
					{
						if(upOrDown == 1)
						convex_changeL.pop_back();
						if(upOrDown == 0)
						convex_changeU.pop_back();
					}

				}
				timeEnd=endIndex+1;
				convex_upper.clear();
				convex_lower.clear();
				startTest = clock();

				double streamdata;
				int sIndexOffset=0;
				long lindex = startIndex;
				Point point1, point2, tpoint1, tpoint2;
				while(lindex<=timeEnd-1)
				{
					sIndexOffset++;
					streamdata=initialStream[lindex-1].y;
					pointCreate(point1,lindex,streamdata+sigma);//upper point
					pointCreate(point2,lindex,streamdata-sigma);//lower point
					if(sIndexOffset == 1)
					{
						if(lindex<initialStream.size())
						{
							lindex++;
							sIndexOffset++;
							streamdata=initialStream[lindex-1].y;
							pointCreate(tpoint1,lindex,streamdata+sigma);//upper point
							pointCreate(tpoint2,lindex,streamdata-sigma);//lower point
							convexInitialization(convex_upper,convex_lower,point1,point2,tpoint1,tpoint2);
							lindex++;
							continue;
						}
					}
					int aa = feasibleRegionUpdate(convex_upper,convex_lower,point1,point2);
					lindex++;
				}

				finishTest=clock();
				lTimeTest = lTimeTest + (double)(finishTest - startTest)*1000;
				finishTest1 = clock();
				ltimeTest1 = ltimeTest1 + (double)(finishTest1 - startTest1)*1000;
				if(upOrDown==0)
					currentUpOrDown=1;
				if(upOrDown==1)
					currentUpOrDown=0;
			}
			//此处包含isLongIntersect==0更新后的，以及isLongIntersect==1的情况
			//判断当前分割的key extremes line是否是上述相交的
			//判断方法：判断方法看当前分割的curentupordown是否和upordown相同，如果相同，存入当前的key extremes line和上一个extremes line的交，如果不同，进行更新操作
			Line linekey;
			Point intersectPoint,updatePointS, updatePointX;
			if(upOrDown==0)
			lineCreate(linekey,convex_lower.back(), convex_upper.front());
			if(upOrDown==1)
			lineCreate(linekey,convex_upper.back(), convex_lower.front());
			intersectPoint=lineIntersection(keyLine,linekey);
			for(int i=startIndex-1; i>=(int)(intersectPoint.x+1);i--)
			{
				//coutDque(convex_upper,convex_lower);
				pointCreate(updatePointS,i,initialStream[i-1].y+sigma);//upper point
				pointCreate(updatePointX,i,initialStream[i-1].y-sigma);//lower point
				int itemp = feasibleRegionUpdateAddForward(convex_upper, convex_lower, updatePointS, updatePointX);
			}

			if(lastseg == true)
			{
				currentUpOrDown=upOrDown;
				blwhile=false;
			}
			if(currentUpOrDown==upOrDown)
			{
				if(currentUpOrDown==1)
				{
					//存入连接点
					Line lineL;
					lineCreate(lineL,convex_upper.back(), convex_lower.front());
					Point pttemp = lineIntersection(keyLine,lineL);
					connectPoint.push_back(lineIntersection(keyLine,lineL));
					numberSeg++;
					if(numberSeg % 257 == 0)
					 {
						 connectFile << endl;
						 connectIndexFile << endl;
					 }
					connectIndexFile <<lineIntersection(keyLine,lineL).x << "    ";
					connectFile << lineIntersection(keyLine,lineL).y << "    ";
					//更新key extremes line，key point
					keyLine=lineL;
					keyPoint=convex_upper.back();
				}
				if(currentUpOrDown==0)
				{
					//存入连接点
					Line lineU;
					lineCreate(lineU,convex_lower.back(), convex_upper.front());
					Point pttemp = lineIntersection(keyLine,lineU);
					connectPoint.push_back(lineIntersection(keyLine,lineU));
					numberSeg++;
					if(numberSeg % 257 == 0)
					 {
						 connectFile << endl;
						 connectIndexFile << endl;
					 }
					 connectIndexFile <<lineIntersection(keyLine,lineU).x << "    ";
					connectFile << lineIntersection(keyLine,lineU).y << "    ";
					//更新key extremes line，key point
					keyLine=lineU;
					keyPoint=convex_lower.back();
				}
				//加上最后的分割的压缩点信息
				if(lastseg)
				{
					connectPoint.push_back(pointCreate(initialStream.size(),keyLine.slope * (initialStream.size()) + keyLine.offset));
					numberSeg++;
					if(numberSeg % 257 == 0)
					 {
						 connectFile << endl;
						 connectIndexFile << endl;
					 }
					 connectIndexFile <<initialStream.size() << "    ";
					connectFile << keyLine.slope * (initialStream.size()) + keyLine.offset << "    ";
				}
				startIndex=timeEnd;
			}
			else
			{
				Point updatePointXI,updatePointSI,intersectPointTemp;
				deque<Point> convex_upper_temp, convex_lower_temp;
				Point maxkeyConvexPoint;
				Line keyLineTemp;
				if(upOrDown==0)
				{
					for(int i=(int)(intersectPoint.x);i>=keyPoint.x;i--)
					{
						lineCreate(linekey,convex_upper.back(), convex_lower.front());
						intersectPointTemp=lineIntersection(keyLine,linekey);
						if(i==(int)(intersectPoint.x))
						{
							if(intersectPointTemp.x >= i && intersectPointTemp.x<(int)(intersectPoint.x+1))
							{
								maxkeyConvexPoint=convex_upper.back();
								keyLineTemp = linekey;
								break;
							}
							else
							{
								pointCreate(updatePointS,i,initialStream[i-1].y+sigma);//upper point
								pointCreate(updatePointX,i,initialStream[i-1].y-sigma);//lower point
								int itemp = feasibleRegionUpdateAddForward(convex_upper, convex_lower, updatePointS, updatePointX);
								if(i==keyPoint.x)
								{
									lineCreate(linekey,convex_upper.back(), convex_lower.front());
									intersectPointTemp=lineIntersection(keyLine,linekey);
									maxkeyConvexPoint=convex_upper.back();
									keyLineTemp = linekey;
									if(intersectPointTemp.x >= i && intersectPointTemp.x < i+1)
									{
										maxkeyConvexPoint=convex_upper.back();
										keyLineTemp = linekey;
										break;
									}
								}
							}
						}
						else
						{
							if(intersectPointTemp.x >= i && intersectPointTemp.x < i+1)
							{
								maxkeyConvexPoint=convex_upper.back();
								keyLineTemp = linekey;
								break;
							}
							else
							{
								pointCreate(updatePointS,i,initialStream[i-1].y+sigma);//upper point
								pointCreate(updatePointX,i,initialStream[i-1].y-sigma);//lower point
								int itemp = feasibleRegionUpdateAddForward(convex_upper, convex_lower, updatePointS, updatePointX);
							}

							if(i==keyPoint.x)
							{
								lineCreate(linekey,convex_upper.back(), convex_lower.front());
								intersectPointTemp=lineIntersection(keyLine,linekey);
								maxkeyConvexPoint=convex_upper.back();
								keyLineTemp = linekey;
								if(intersectPointTemp.x >= i && intersectPointTemp.x < i+1)
								{
									maxkeyConvexPoint=convex_upper.back();
									keyLineTemp = linekey;
									break;
								}
							}
						}

					}
					//压入连接点
					connectPoint.push_back(lineIntersection(keyLineTemp, keyLine));
					numberSeg++;
					if(numberSeg % 257 == 0)
					 {
						 connectFile << endl;
						 connectIndexFile << endl;
					 }
					 connectIndexFile <<lineIntersection(keyLineTemp, keyLine).x << "    ";
					connectFile << lineIntersection(keyLineTemp, keyLine).y << "    ";
					//更新
					keyLine=keyLineTemp;
					keyPoint=maxkeyConvexPoint;
					startIndex=timeEnd;
				}

				if(upOrDown==1)
				{
					for(int i=(int)(intersectPoint.x);i>=keyPoint.x;i--)
					{
						lineCreate(linekey,convex_lower.back(), convex_upper.front());
						intersectPointTemp=lineIntersection(keyLine,linekey);
						if(i==(int)(intersectPoint.x))
						{
							if(intersectPointTemp.x >= i && intersectPointTemp.x<(int)(intersectPoint.x+1))
							{
								maxkeyConvexPoint=convex_lower.back();
								keyLineTemp = linekey;
								break;
							}
							else
							{
								pointCreate(updatePointS,i,initialStream[i-1].y+sigma);//upper point
								pointCreate(updatePointX,i,initialStream[i-1].y-sigma);//lower point
								int itemp = feasibleRegionUpdateAddForward(convex_upper, convex_lower, updatePointS, updatePointX);
							}
							if(i==keyPoint.x)
							{
								lineCreate(linekey,convex_lower.back(), convex_upper.front());
								intersectPointTemp=lineIntersection(keyLine,linekey);
								maxkeyConvexPoint=convex_lower.back();
								keyLineTemp = linekey;
								if(intersectPointTemp.x >= i && intersectPointTemp.x < i+1)
								{
									maxkeyConvexPoint=convex_lower.back();
									keyLineTemp = linekey;
									break;
								}
							}
						}
						else
						{
							if(intersectPointTemp.x >= i && intersectPointTemp.x < i+1)
							{
								maxkeyConvexPoint=convex_lower.back();
								keyLineTemp = linekey;
								break;
							}
							else
							{
								pointCreate(updatePointS,i,initialStream[i-1].y+sigma);//upper point
								pointCreate(updatePointX,i,initialStream[i-1].y-sigma);//lower point
								int itemp = feasibleRegionUpdateAddForward(convex_upper, convex_lower, updatePointS, updatePointX);
							}

							if(i==keyPoint.x)
							{
								lineCreate(linekey,convex_lower.back(), convex_upper.front());
								intersectPointTemp=lineIntersection(keyLine,linekey);
								maxkeyConvexPoint=convex_lower.back();
								keyLineTemp = linekey;
								if(intersectPointTemp.x >= i && intersectPointTemp.x < i+1)
								{
									maxkeyConvexPoint=convex_lower.back();
									keyLineTemp = linekey;
									break;
								}
							}
						}
					}
					//压入连接点
					connectPoint.push_back(lineIntersection(keyLineTemp, keyLine));
					numberSeg++;
					if(numberSeg % 257 == 0)
					{
						connectFile << endl;
						connectIndexFile << endl;
					}
					connectIndexFile <<lineIntersection(keyLineTemp, keyLine).x << "    ";
					connectFile << lineIntersection(keyLineTemp, keyLine).y << "    ";
					//更新
					keyLine=keyLineTemp;
					keyPoint=maxkeyConvexPoint;
					startIndex=timeEnd;
				}
			}
		}
		upOrDown=currentUpOrDown;
		convex_lower.clear();
		convex_upper.clear();
	}
	finish = clock();
	duration = (double)(finish - start)*1000;
	resultFile<<"time cost: "<< duration / 1000 <<" ms"<<endl;
	resultFile<<"The number of segment is: "<<connectPoint.size()<<endl;

	//计算总体误差
	double sumError=0;
	deque<Point>::iterator iter = connectPoint.begin();
	for (;iter != connectPoint.end(); iter ++)
	{
		Point pBegin = (*iter);
		iter++;
		Point pEnd = (*iter);
		if(pBegin.x == 1)
			for(long i=pBegin.x; i<=pEnd.x;i++)
			{
				Line l = lineCreate(pBegin, pEnd);
				double iValue = l.slope * i + l.offset;
				double dPlus = iValue - initialStream[i-1].y;
				sumError += dPlus * dPlus;
			}
		else
			for(long i=pBegin.x+1; i<=pEnd.x;i++)
			{
				Line l = lineCreate(pBegin, pEnd);
				double iValue = l.slope * i + l.offset;
				double dPlus = iValue - initialStream[i-1].y;
				sumError += dPlus * dPlus;
			}
		iter--;
	}
	sumError = sqrt(sumError/initialStream.size());
	resultFile<<"sumError: "<<sumError<<endl;

	streamFile.close();
	resultFile.close();
	connectIndexFile.close();
	connectFile.close();
	return;
}

void optimalProcessingFS(string fileName, string resultName, string connectionIndex, string connectionValue, double sigma)
{
	ifstream streamFile(fileName.c_str());
	ofstream resultFile(resultName.c_str(),ios::app);
	ofstream connectIndexFile(connectionIndex.c_str());
	ofstream connectFile(connectionValue.c_str());

	//装入所有数据
	deque<Point> initialStream;
	Point initialPoint;
	initialPoint=pointCreate(0,0);
	long initialIndex=1;
	double initialData;
	double maxVal, minVal;
	maxVal = -10000000000;
	minVal = 10000000000;

	FILE *fp;
	fp=fopen(fileName.c_str(),"r");
	if (fp == NULL)
	{
	    exit(1);
	}
	while(!feof(fp))
	{
	     fscanf(fp,"%lf",&initialData);
	     pointCreate(initialPoint,initialIndex,initialData);
	     if(maxVal < initialData)
	     {
	     	maxVal = initialData;
	     }
	     if(minVal > initialData)
	     {
	     	minVal = initialData;
	     }
	     initialStream.push_back(initialPoint);
	     initialIndex++;
	}
	fclose(fp);

	resultFile<<"误差："<< " "<< sigma << endl;
	resultFile<<""<< " " << endl;
	deque<Point> convex_upper, convex_lower; //convex_upper上极限点的容器
	long sIndex; //正在处理第几个数据
	long startIndex; //分割的第一个点
	long sIndexOffset; //分割内的第几个数据
	long numberSeg; //分割次数
	Point point1, point2, tpoint1, tpoint2;
	Line tempU, minLast, maxLast;
	double streamdata;
	deque<Point> connectPoint;

	//初始化所有用到的数据
	tempU=lineCreate(0,0);
	streamdata=0;
	sIndex = 0;
	startIndex = 0;
	sIndexOffset = 0;
	numberSeg = 0;

	//时间消耗
	clock_t start,finish;
	double duration;
	start = clock();
	Point p1, p2,startPoint;
	Line l1, l2;
	int startIndexs = 1;
	FS newFS;
	FS lastFS;
	bool b=true;
	//开始处理
	while(sIndex<initialStream.size())
	{
		sIndex++;
		if(b==true)
		{
			if(numberSeg == 0)
			{
				//初始话可行空间
				startPoint = initialStream[sIndex-1];
				connectPoint.push_back(startPoint);
				connectIndexFile <<startPoint.x << "    ";
				connectFile << startPoint.y << "    ";
				sIndex++;
				lastFS.lowerBound = (initialStream[sIndex - 1].y - startPoint.y - sigma)
						/ (initialStream[sIndex - 1].x - startPoint.x);
				lastFS.upperBound = (initialStream[sIndex - 1].y - startPoint.y + sigma)
								/ (initialStream[sIndex - 1].x - startPoint.x);
				b = false;
				continue;
			}
		}
		else
		{
			//更新FS
			newFS.lowerBound = (initialStream[sIndex - 1].y - startPoint.y - sigma)
							/ (initialStream[sIndex - 1].x - startPoint.x);
			newFS.upperBound = (initialStream[sIndex - 1].y - startPoint.y + sigma)
							/ (initialStream[sIndex - 1].x - startPoint.x);

			if (newFS.lowerBound >= lastFS.upperBound || newFS.upperBound <= lastFS.lowerBound) {
				Line lineTemp1 = lineCreate(lastFS.upperBound,startPoint);
				Line lineTemp2 = lineCreate(lastFS.lowerBound,startPoint);
				numberSeg++;
				double p = 0;
				p = (lineTemp1.slope + lineTemp2.slope)/2;
				startPoint = pointCreate(sIndex - 1, p * (sIndex - 1) + startPoint.y - p * startPoint.x);
				if(numberSeg % 257 == 0)
				{
					connectFile << endl;
					connectIndexFile << endl;
				}
				connectIndexFile <<startPoint.x << "    ";
				connectFile << startPoint.y << "    ";
				connectPoint.push_back(startPoint);
				lastFS.lowerBound = (initialStream[sIndex - 1].y - startPoint.y - sigma)
						/ (initialStream[sIndex - 1].x - startPoint.x);
				lastFS.upperBound = (initialStream[sIndex - 1].y - startPoint.y + sigma)
								/ (initialStream[sIndex - 1].x - startPoint.x);
				continue;
			}
			if (newFS.lowerBound > lastFS.lowerBound) { // narrow fs
				lastFS.lowerBound = newFS.lowerBound;
			}
			if (newFS.upperBound < lastFS.upperBound) {
				lastFS.upperBound = newFS.upperBound;
			}
			continue;

		}
	}
	if(sIndex>0)
	{
		Line lineTemp1 = lineCreate(lastFS.upperBound,startPoint);
		Line lineTemp2 = lineCreate(lastFS.lowerBound,startPoint);
		numberSeg++;
		double p = 0;
		p = (lineTemp1.slope + lineTemp2.slope)/2;
		startPoint = pointCreate(sIndex - 1, p * (sIndex - 1) + startPoint.y - p * startPoint.x);
		if(numberSeg % 257 == 0)
		{
			connectFile << endl;
			connectIndexFile << endl;
		}
		connectIndexFile <<startPoint.x << "    ";
		connectFile << startPoint.y << "    ";
		connectPoint.push_back(startPoint);
	}

	finish = clock();

	duration = (double)(finish - start)*1000;
	resultFile<<"time cost: "<<duration / 1000<<" ms"<<endl;
	resultFile<<"The number of segment is: "<<numberSeg<<endl;

	//计算总体误差
	double sumError=0;
	deque<Point>::iterator iter = connectPoint.begin();
	for (;iter != connectPoint.end(); iter ++)
	{
		Point pBegin = (*iter);
		iter++;
		Point pEnd = (*iter);
		if(pBegin.x == 1)
			for(long i=pBegin.x; i<=pEnd.x;i++)
			{
				Line l = lineCreate(pBegin, pEnd);
				double iValue = l.slope * i + l.offset;
				double dPlus = iValue - initialStream[i-1].y;
				sumError += dPlus * dPlus;
			}
		else
			for(long i=pBegin.x+1; i<=pEnd.x;i++)
			{
				Line l = lineCreate(pBegin, pEnd);
				double iValue = l.slope * i + l.offset;
				double dPlus = iValue - initialStream[i-1].y;
				sumError += dPlus * dPlus;
			}
		iter--;
	}
	sumError = sqrt(sumError / initialStream.size());
	resultFile<<"sumErrorFSW: "<<sumError<<endl;


	streamFile.close();
	resultFile.close();
	connectIndexFile.close();
	connectFile.close();
	return;
}

