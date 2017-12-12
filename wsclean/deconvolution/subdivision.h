#ifndef SUBDIVISION_H
#define SUBDIVISION_H

#include <algorithm>
#include <set>

class Subdivision
{
public:
	Subdivision(size_t width, size_t height) :
		_width(width), _height(height)
	{
	}
	
	struct Coord {
		Coord() = default;
		Coord(size_t _x, size_t _y) : x(_x), y(_y) {}
		size_t x, y;
	};
	
	struct Visit {
		double distance;
		Coord to, from;
		bool operator<(const Visit& rhs) const {
			return distance < rhs.distance;
		}
	};
	
	void DivideVertically(const double* image, double* output, size_t x1, size_t x2)
	{
		using visitset = std::multiset<Visit>;
		visitset visits;
		
		for(size_t x=x1; x!=x2; ++x)
		{
			Visit visit;
			visit.distance = 0.0;
			visit.to = Coord(x, 0);
			visit.from = Coord(x, 0);
			visits.insert(visit);
		}
		ao::uvector<Coord> path((x2-x1) * _height);
		std::fill(output, output+_width*_height, std::numeric_limits<double>::max());
		Visit visit;
		while(!visits.empty())
		{
			visit = *visits.begin();
			visits.erase(visits.begin());
			size_t x = visit.to.x, y = visit.to.y;
			double curDistance = output[x + y*_width];
			double newDistance = visit.distance + std::fabs(image[x + y*_width]);
			//std::cout << x << ',' << y << " " << curDistance << " " << newDistance << '\n';
			if(newDistance < curDistance)
			{
				output[x + y*_width] = newDistance;
				path[x-x1 + y*(x2-x1)] = visit.from;
				if(y == _height-1)
					break;
				visit.distance = newDistance;
				visit.from = visit.to;
				if(x > x1) {
					visit.to = Coord(x-1, y+1);
					visits.insert(visit);
					visit.to = Coord(x-1, y);
					visits.insert(visit);
				}
				visit.to = Coord(x, y+1);
				visits.insert(visit);
				if(x < x2-1) {
					visit.to = Coord(x+1, y+1);
					visits.insert(visit);
					visit.to = Coord(x+1, y);
					visits.insert(visit);
				}
			}
		}
		std::fill(output, output+_width*_height, 0.0);
		Coord pCoord = visit.to;
		while(pCoord.y > 0) {
			output[pCoord.x + pCoord.y*_width] = 1.0;
			pCoord = path[pCoord.x-x1 + pCoord.y*(x2-x1)];
		}
	}
	
	void DivideHorizontally(const double* image, double* output, size_t y1, size_t y2)
	{
		using visitset = std::multiset<Visit>;
		visitset visits;
		
		for(size_t y=y1; y!=y2; ++y)
		{
			Visit visit;
			visit.distance = 0.0;
			visit.to = Coord(0, y);
			visit.from = Coord(0, y);
			visits.insert(visit);
		}
		ao::uvector<Coord> path(_width * (y2-y1));
		std::fill(output, output+_width*_height, std::numeric_limits<double>::max());
		Visit visit;
		while(!visits.empty())
		{
			visit = *visits.begin();
			visits.erase(visits.begin());
			size_t x = visit.to.x, y = visit.to.y;
			double curDistance = output[x + y*_width];
			double newDistance = visit.distance + std::fabs(image[x + y*_width]);
			//std::cout << x << ',' << y << " " << curDistance << " " << newDistance << '\n';
			if(newDistance < curDistance)
			{
				output[x + y*_width] = newDistance;
				path[x + (y-y1)*_width] = visit.from;
				if(x == _width-1)
					break;
				visit.distance = newDistance;
				visit.from = visit.to;
				if(y > y1) {
					visit.to = Coord(x+1, y-1);
					visits.insert(visit);
					visit.to = Coord(x, y-1);
					visits.insert(visit);
				}
				visit.to = Coord(x+1, y);
				visits.insert(visit);
				if(y < y2-1) {
					visit.to = Coord(x+1, y+1);
					visits.insert(visit);
					visit.to = Coord(x, y+1);
					visits.insert(visit);
				}
			}
		}
		std::fill(output, output+_width*_height, 0.0);
		Coord pCoord = visit.to;
		while(pCoord.x > 0) {
			output[pCoord.x + pCoord.y*_width] = 1.0;
			pCoord = path[pCoord.x + (pCoord.y-y1)*_width];
		}
	}
	
	struct HorizontalSplits {
		size_t leftX, leftY;
		size_t rightX, rightY;
		size_t minLeftY, maxLeftY;
		size_t minRightY, maxRightY;
	};
	
	HorizontalSplits DivideHorizontallyWithSplits(const double* image, const double* verticalBoundary, double* output, size_t y1, size_t y2)
	{
		using visitset = std::multiset<Visit>;
		visitset visits;
		
		for(size_t y=y1; y!=y2; ++y)
		{
			Visit visit;
			visit.distance = 0.0;
			visit.to = Coord(0, y);
			visit.from = Coord(0, y);
			visits.insert(visit);
		}
		ao::uvector<Coord> path(_width * (y2-y1));
		std::fill(output, output+_width*_height, std::numeric_limits<double>::max());
		Visit visit;
		while(!visits.empty())
		{
			visit = *visits.begin();
			visits.erase(visits.begin());
			size_t x = visit.to.x, y = visit.to.y;
			double curDistance = output[x + y*_width];
			double newDistance = visit.distance + std::fabs(image[x + y*_width]);
			//std::cout << x << ',' << y << " " << curDistance << " " << newDistance << '\n';
			if(newDistance < curDistance)
			{
				output[x + y*_width] = newDistance;
				path[x + (y-y1)*_width] = visit.from;
				if(x == _width-1)
					break;
				visit.distance = newDistance;
				visit.from = visit.to;
				if(y > y1) {
					visit.to = Coord(x+1, y-1);
					visits.insert(visit);
					visit.to = Coord(x, y-1);
					visits.insert(visit);
				}
				visit.to = Coord(x+1, y);
				visits.insert(visit);
				if(y < y2-1) {
					visit.to = Coord(x+1, y+1);
					visits.insert(visit);
					visit.to = Coord(x, y+1);
					visits.insert(visit);
				}
			}
		}
		std::fill(output, output+_width*_height, 0.0);
		Coord pCoord = visit.to;
		HorizontalSplits splits;
		splits.minLeftY = y2;
		splits.maxLeftY = y1;
		splits.minRightY = y2;
		splits.maxRightY = y1;
		bool hasRight = false;
		while(pCoord.x > 0)
		{
			output[pCoord.x + pCoord.y*_width] = 1.0;
			if(!hasRight)
			{
				splits.minRightY = std::min(splits.minRightY, pCoord.y);
				splits.maxRightY = std::max(splits.maxRightY, pCoord.y);
			}
			if(verticalBoundary[pCoord.x + pCoord.y*_width] != 0.0)
			{
				splits.leftX = pCoord.x;
				splits.leftY = pCoord.y;
				if(!hasRight)
				{
					splits.rightX = pCoord.x;
					splits.rightY = pCoord.y;
					hasRight = true;
				}
				// Sometimes it happens that the boundaries take different paths --
				// if so make sure to reset the left boundaries again.
				splits.minLeftY = y2;
				splits.maxLeftY = y1;
			}
			else {
				if(hasRight)
				{
					splits.minLeftY = std::min(splits.minLeftY, pCoord.y);
					splits.maxLeftY = std::max(splits.maxLeftY, pCoord.y);
				}
			}
			pCoord = path[pCoord.x + (pCoord.y-y1)*_width];
		}
		splits.minLeftY = std::min(splits.leftY, splits.minLeftY);
		splits.maxLeftY = std::max(splits.leftY, splits.maxLeftY);
		splits.minRightY = std::min(splits.rightY, splits.minRightY);
		splits.maxRightY = std::max(splits.rightY, splits.maxRightY);
		return splits;
	}
private:
	size_t _width, _height;
};

#endif

