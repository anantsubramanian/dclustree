#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <unistd.h>

using namespace std;

#define INITIALCLUSTERS 7
#define INITIALPOINTS 30
#define NUMCLUSTERS 25
//#define NUMPOINTS 1000000
#define LOWERX -200
#define UPPERX 200
#define LOWERY -200
#define UPPERY 200
#define MINRADIUS 5
#define RADIUS 20
#define MINSPEED 25000
#define MAXSPEED 250000
#define LOWERPOINTS 75
#define UPPERPOINTS 200

class Location
{
	public:
		int x, y;
		Location(int x, int y)
		{
			this->x = x;
			this->y = y;
		}
};

vector<Location> locations;

int main(int argc, char *argv[])
{
	int curx, cury, radius;
	int usleepamount;
	srand(time(NULL));
	for(int i = 0; i < NUMCLUSTERS; i++)
	{
		curx = LOWERX + rand()%(UPPERX - LOWERX + 1);
		cury = LOWERY + rand()%(UPPERY - LOWERY + 1);
		locations.push_back(Location(curx, cury));
	}
	usleepamount = MINSPEED + rand()%(MAXSPEED - MINSPEED + 1);
	int changeafter = LOWERPOINTS + rand()%(UPPERPOINTS - LOWERPOINTS + 1);
	int count = 0;
	int initialized = 0;
	int numpoints;
	cin>>numpoints;
	for(int i = 0; i < numpoints; i++)
	{
		int randomcluster = rand()%((int)locations.size());
		if (initialized < INITIALCLUSTERS)
		{
			for(int j = 0; j < INITIALPOINTS; j++)
			{
				int radiusx = MINRADIUS + rand()%(RADIUS - MINRADIUS + 1);
				int radiusy = MINRADIUS + rand()%(RADIUS - MINRADIUS + 1);
				int deltax = rand()%RADIUS;
				int x = locations[randomcluster].x + (rand()%2 == 0 ? radiusx : -radiusx);
				int y = locations[randomcluster].y + (rand()%2 == 0 ? radiusy : -radiusy);
				printf("%d %d\n", x, y);
			}
			initialized++;
			i--;
		}
		else
		{
			int radiusx = MINRADIUS + rand()%(RADIUS - MINRADIUS + 1);
			int radiusy = MINRADIUS + rand()%(RADIUS - MINRADIUS + 1);
			int deltax = rand()%RADIUS;
			int x = locations[randomcluster].x + (rand()%2 == 0 ? radiusx : -radiusx);
			int y = locations[randomcluster].y + (rand()%2 == 0 ? radiusy : -radiusy);
			printf("%d %d\n", x, y);
		}
		//usleep(usleepamount);
		count++;
		if (count == changeafter)
		{
			count = 0;
			usleepamount = MINSPEED + rand()%(MAXSPEED - MINSPEED + 1);
			changeafter = LOWERPOINTS + rand()%(UPPERPOINTS - LOWERPOINTS + 1);
		}
	}
	return 0;		
}
