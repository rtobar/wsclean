#ifndef SUBDIVISION_H
#define SUBDIVISION_H

#include <algorithm>
#include <queue>

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
			return distance > rhs.distance;
		}
	};
	
	void DivideVertically(const double* image, double* output, size_t x1, size_t x2) const
	{
		using visitset = std::priority_queue<Visit>;
		visitset visits;
		
		for(size_t x=x1; x!=x2; ++x)
		{
			Visit visit;
			visit.distance = 0.0;
			visit.to = Coord(x, 0);
			visit.from = Coord(x, 0);
			visits.push(visit);
		}
		ao::uvector<Coord> path((x2-x1) * _height);
		std::fill(output, output+_width*_height, std::numeric_limits<double>::max());
		Visit visit;
		while(!visits.empty())
		{
			visit = visits.top();
			visits.pop();
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
					visits.push(visit);
					visit.to = Coord(x-1, y);
					visits.push(visit);
				}
				visit.to = Coord(x, y+1);
				visits.push(visit);
				if(x < x2-1) {
					visit.to = Coord(x+1, y+1);
					visits.push(visit);
					visit.to = Coord(x+1, y);
					visits.push(visit);
				}
			}
		}
		std::fill(output, output+_width*_height, 0.0);
		Coord pCoord = visit.to;
		while(pCoord.y > 0) {
			output[pCoord.x + pCoord.y*_width] = 1.0;
			pCoord = path[pCoord.x-x1 + pCoord.y*(x2-x1)];
		}
		output[pCoord.x] = 1.0;
	}
	
	void DivideHorizontally(const double* image, double* output, size_t y1, size_t y2) const
	{
		using visitset = std::priority_queue<Visit>;
		visitset visits;
		
		for(size_t y=y1; y!=y2; ++y)
		{
			Visit visit;
			visit.distance = 0.0;
			visit.to = Coord(0, y);
			visit.from = Coord(0, y);
			visits.push(visit);
		}
		ao::uvector<Coord> path(_width * (y2-y1));
		std::fill(output, output+_width*_height, std::numeric_limits<double>::max());
		Visit visit;
		while(!visits.empty())
		{
			visit = visits.top();
			visits.pop();
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
					visits.push(visit);
					visit.to = Coord(x, y-1);
					visits.push(visit);
				}
				visit.to = Coord(x+1, y);
				visits.push(visit);
				if(y < y2-1) {
					visit.to = Coord(x+1, y+1);
					visits.push(visit);
					visit.to = Coord(x, y+1);
					visits.push(visit);
				}
			}
		}
		std::fill(output, output+_width*_height, 0.0);
		Coord pCoord = visit.to;
		while(pCoord.x > 0) {
			output[pCoord.x + pCoord.y*_width] = 1.0;
			pCoord = path[pCoord.x + (pCoord.y-y1)*_width];
		}
		output[pCoord.y*_width] = 1.0;
	}
	
	struct HorizontalSplits {
		size_t leftX, leftY;
		size_t rightX, rightY;
		size_t minLeftY, maxLeftY;
		size_t minRightY, maxRightY;
	};
	
	HorizontalSplits DivideHorizontallyWithSplits(const double* image, const double* verticalBoundary, double* output, size_t y1, size_t y2) const
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
	
	void GetBoundingMask(const double* subdivision, size_t subImgX, size_t subImgY,
		bool* mask, bool* visited, size_t& x, size_t& y, size_t& subWidth, size_t& subHeight) const
	{
		// This function performs a flood fill of the subimage. Two results are stored: the
		// mask containing all the pixels that belong to this subimage, as well as
		// a visited image that contains the union of all masks so far created.
		// The latter is necessary to prevent storing edge pixels in two subimages.
		// (While this could also be deduced by tracing the subdivision paths, this
		// was simpler -- at the cost of one extra array of booleans)
		
		std::fill(mask, mask+_width*_height, false);
		std::vector<std::pair<size_t,size_t>> stack;
		stack.emplace_back(subImgX, subImgY);
		x = _width;
		y = _height;
		size_t x2 = 0, y2 = 0;
		while(!stack.empty())
		{
			size_t
				curX = stack.back().first,
				curY = stack.back().second;
			stack.pop_back();
			if(!visited[curX + curY*_width])
			{
				visited[curX + curY*_width] = true;
				mask[curX + curY*_width] = true;
				x = std::min(x, curX); x2 = std::max(x2, curX);
				y = std::min(y, curY); y2 = std::max(y2, curY);
				if(subdivision[curX + curY * _width] == 0.0)
				{
					if(curX > 0)
						stack.emplace_back(curX-1, curY);
					// When moving down / to the right, we do not allow to be on the border,
					// because the right border is excluded from the subdivision edge.
					// (note that the 'visited' mask prevents overlap, but this way
					// subimages consistenly include the left edge and exclude the right edge
					// as long as the edge is not surrounded on three sides, in which case it is
					// more or less random)
					if(curX < _width-1 && subdivision[curX + 1 + curY * _width] == 0.0)
						stack.emplace_back(curX+1, curY);
					if(curY > 0)
						stack.emplace_back(curX, curY-1);
					if(curY < _height-1 && subdivision[curX + (curY+1) * _width] == 0.0)
						stack.emplace_back(curX, curY+1);
				}
			}
		}
		if(x2 < x)
		{
			subWidth = 0;
			subHeight = 0;
		}
		else {
			subWidth = x2 - x;
			subHeight = y2 - y;
		}
		// If dimensions start of even, keep subimages even too
		if(_width%2 == 0)
		{
			if(subWidth%2 != 0)
			{
				if(subWidth + x >= _width)
					--x;
				++subWidth;
			}
		}
		if(_height%2 == 0)
		{
			if(subHeight%2 != 0)
			{
				if(subHeight + y >= _height)
					--y;
				++subHeight;
			}
		}
	}
private:
	size_t _width, _height;
};

#endif

