#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>
#include <mpi.h>

using namespace std;

#define LAMBDA 0.00001
#define BETA 2
#define NUMINSERTS 300
#define MASTER 0
#define SLAVE1 1
#define SLAVE2 2
#define INITIALSIZE 10
#define UPDATECHECK 30

int points = 0;

class Point
{
	public:
		int timestamp;
		int x, y;
		Point(int x, int y, int timestamp)
		{
			this->x = x;
			this->y = y;
			this->timestamp = timestamp;
		}

};

class CF
{
	public:
		double n;
		double lsx, ssx, lsy, ssy;
		int timestamp;
		
		CF(int lsx, int lsy, int timestamp, double n)
		{
			this->n = n;
			this->lsx = lsx;
			this->lsy = lsy;
			this->ssx = this->lsx * this->lsx;
			this->ssy = this->lsy * this->lsy;
			this->timestamp = timestamp;
		}

		CF(double n, double lsx, double lsy, double ssx, double ssy, int timestamp)
		{
			this->n = n;
			this->lsx = lsx;
			this->lsy = lsy;
			this->ssx = ssx;
			this->ssy = ssy;
			this->timestamp = timestamp;
		}

		void update(int timestamp)
		{
			this->n *= pow(BETA,-LAMBDA*(timestamp - this->timestamp));
			this->lsx *= pow(BETA,-LAMBDA*(timestamp - this->timestamp));
			this->lsy *= pow(BETA,-LAMBDA*(timestamp - this->timestamp));
			this->ssx *= pow(BETA,-LAMBDA*(timestamp - this->timestamp));
			this->ssy *= pow(BETA,-LAMBDA*(timestamp - this->timestamp));
			this->timestamp = timestamp;
		}

};

class Node
{
	public:
		bool isleaf;
		int size;
		CF **cf;
		Node **child;
		Node *parent;
		
		Node(Node *parent, int m, int M, int l, int L, bool isleaf)
		{
			this->size = 0;
			this->cf = new CF*[M+1];
			this->child = new Node*[M+1];
			for(int i = 0; i < M+1; i++)
			{
				this->cf[i] = NULL;
				this->child[i] = NULL;
			}
			this->parent = parent;
			this->isleaf = isleaf;
		}

		void addCF(CF *cf, Node *relatedchild)
		{
			// Invoked during splits
			this->cf[this->size] = cf;
			this->child[this->size] = relatedchild;
			this->size++;
		}

		void addCF(double n, double lsx, double lsy, double ssx, double ssy, int timestamp, Node *relatedchild)
		{
			this->cf[this->size] = new CF(n, lsx, lsy, ssx, ssy, timestamp);
			this->child[this->size] = relatedchild;
			this->size++;
		}

		void addPoint(Point *p)
		{
			// Only called on a leaf node instance
			// child pointers remain NULL
			this->cf[this->size] = new CF(p->x, p->y, p->timestamp, 1);
			this->size++;
		}

		CF* getCF(int newtimestamp)
		{
			double n = 0.0, lsx = 0.0, lsy = 0.0, ssx = 0.0, ssy = 0.0;
			for(int i = 0; i < this->size; i++)
			{
				n += this->cf[i]->n;
				lsx += this->cf[i]->lsx;
				lsy += this->cf[i]->lsy;
				ssx += this->cf[i]->ssx;
				ssy += this->cf[i]->ssy;
			}

			return new CF(n, lsx, lsy, ssx, ssy, newtimestamp);
		}

};

class ClusTree
{
	int m, M, l, L;
	Node *root;
	int deltaT, lastpointtime;
	public:
		ClusTree(int m, int M, int l, int L)
		{
			this->m = m;
			this->M = M;
			this->l = l;
			this->L = L;
			this->root = NULL;
			this->deltaT = 0;
			this->lastpointtime = 0;
		}

		double getDistance(Point *p, CF *cf)
		{
			// Get Euclidean distance of point from centroid of CF
			return sqrt(pow((p->x - (cf->lsx/cf->n)), 2) + pow((p->y - (cf->lsy/cf->n)), 2));
		}

		double getDistance(CF *cf1, CF *cf2)
		{
			double x1, y1, x2, y2;
			x1 = cf1->lsx / cf1->n;
			x2 = cf2->lsx / cf2->n;
			y1 = cf1->lsy / cf1->n;
			y2 = cf2->lsy / cf2->n;
			return sqrt(pow(x1-x2, 2) + pow(y1-y2,2));
		}
		
		bool needsSplit(Node *node)
		{
			if (node->size <= M)
				return false;
			return true;
		}
		
		bool checkedall(vector<int> &positions)
		{
			for(int i = (M+1)/2 - 1, j = M; i >=	0; i--, j--)
				if(positions[i] != j)
					return false;
			return true;
		}
		
		double getIntraGroupDistance(Node *node, vector<int> &group1, vector<int> &group2)
		{
			double distance = 0.0;
			for(int i = 0; i < group1.size(); i++)
				for(int j = i+1; j < group1.size(); j++)
					distance += getDistance(node->cf[group1[i]], node->cf[group1[j]]);
			for(int i = 0; i < group2.size(); i++)
				for(int j = i+1; j < group2.size(); j++)
					distance += getDistance(node->cf[group2[i]], node->cf[group2[j]]);
			return distance;
		}

		void getNextPosition(vector<int> &positions)
		{
			int i, j;
			for(i = positions.size()-1, j = M; i >= 0; i--, j--)
			{
				if(positions[i] < j)
					break;
			}
			positions[i]++;
			for(int k = i+1; k < positions.size(); k++)
				positions[k] = positions[k-1]+1;
		}

		void split(Node *node, int parentCFpos, int newtimestamp)
		{
			//cout<<"Splitting "<<node->getCF()->lsx<<" "<<node->getCF()->lsy<<"\n";
			vector<int> bestgroup1, bestgroup2;
			vector<int> group1, group2;
			vector<int> positions((M+1)/2);
			for(int i = 0; i < (M+1)/2; i++)
				positions[i] = i;
			group1 = positions;
			for(int i = 0; i < M+1; i++)
				if(find(positions.begin(), positions.end(), i) == positions.end())
					group2.push_back(i);
			bestgroup1 = group1;
			bestgroup2 = group2;
			double mindist = getIntraGroupDistance(node, group1, group2);
			do
			{
				getNextPosition(positions);
				group1 = positions;
				group2 = vector<int>();
				for(int i = 0; i < M+1; i++)
					if(find(positions.begin(), positions.end(), i) == positions.end())
						group2.push_back(i);
				double tempdist = getIntraGroupDistance(node, group1, group2);
				if (tempdist < mindist)
				{
					bestgroup1 = group1;
					bestgroup2 = group2;
					mindist = tempdist;
				}
			} while(!checkedall(positions));
			
			// Create new nodes to accomodate the two groups of cluster features
			Node *n1 = NULL, *n2 = NULL;

			if (node->parent == NULL)
			{
				// At root node, so split should create new root node
				this->root = new Node(NULL, m, M, l, L, false);
				n1 = new Node(this->root, m, M, l, L, node->isleaf);
				n2 = new Node(this->root, m, M, l, L, node->isleaf);
				for(int i = 0; i < bestgroup1.size(); i++)
				{
					n1->addCF(node->cf[bestgroup1[i]], node->child[bestgroup1[i]]);
					if (node->child[bestgroup1[i]])
						node->child[bestgroup1[i]]->parent = n1;
				}
				for(int i = 0; i < bestgroup2.size(); i++)
				{
					n2->addCF(node->cf[bestgroup2[i]], node->child[bestgroup2[i]]);
					if (node->child[bestgroup2[i]])
						node->child[bestgroup2[i]]->parent = n2;
				}
				this->root->addCF(n1->getCF(newtimestamp), n1);
				this->root->addCF(n2->getCF(newtimestamp), n2);
				delete node;
			}
			else
			{
				// Node could be internal or leaf node
				n1 = new Node(node->parent, m, M, l, L, node->isleaf);
				n2 = new Node(node->parent, m, M, l, L, node->isleaf);
				for(int i = 0; i < bestgroup1.size(); i++)
				{
					n1->addCF(node->cf[bestgroup1[i]], node->child[bestgroup1[i]]);
					if (node->child[bestgroup1[i]])
						node->child[bestgroup1[i]]->parent = n1;
				}
				for(int i = 0; i < bestgroup2.size(); i++)
				{
					n2->addCF(node->cf[bestgroup2[i]], node->child[bestgroup2[i]]);
					if (node->child[bestgroup2[i]])
						node->child[bestgroup2[i]]->parent = n2;
				}
				CF *cftodelete = node->parent->cf[parentCFpos];
				node->parent->cf[parentCFpos] = n1->getCF(newtimestamp);
				node->parent->child[parentCFpos] = n1;
				node->parent->addCF(n2->getCF(newtimestamp), n2);
				delete cftodelete;
				delete node;
			}
		}
		
		CF* insert(Point *p, Node *curnode, int whichCFofParent, int newtimestamp)
		{
			if (curnode->isleaf)
			{
				double minimum = curnode->cf[0]->n;
				int minpos = 0;
				for(int i = 1; i < curnode->size; i++)
				{	
					curnode->cf[i]->update(newtimestamp);
					if (curnode->cf[i]->n < minimum)
					{
						minimum = curnode->cf[i]->n;
						minpos = i;
					}
				}
				if (minimum < pow(BETA, -LAMBDA * this->deltaT * NUMINSERTS) && curnode->size+1 > M)
				{
					//cout<<minimum<<" is less than "<<pow(BETA, -LAMBDA * this->deltaT * NUMINSERTS)<<"\n";
					CF *tempcf = curnode->cf[minpos];
					curnode->cf[minpos] = new CF(p->x, p->y, newtimestamp, 1);
					return tempcf;
				}
				else
				{
					curnode->addPoint(p);
					if (needsSplit(curnode))
						split(curnode, whichCFofParent, newtimestamp);
					return NULL;
				}
			}
			else
			{
				// Find closest cluster feature to recurse on while not a child
				for(int i = 0; i < curnode->size; i++)
					curnode->cf[i]->update(newtimestamp);
				double mindist = getDistance(p, curnode->cf[0]);
				int insertpos = 0;
				for (int i = 1; i < curnode->size; i++)
				{
					double tempdist = getDistance(p, curnode->cf[i]);
					if (tempdist < mindist)
					{
						mindist = tempdist;
						insertpos = i;
					}
				}
				// Update CFs on the path to the leaf
				// TODO: Update timestamp
				curnode->cf[insertpos]->lsx += p->x;
				curnode->cf[insertpos]->lsy += p->y;
				curnode->cf[insertpos]->ssx += p->x*p->x;
				curnode->cf[insertpos]->ssy += p->y*p->y;
				curnode->cf[insertpos]->n += 1;

				CF *tempCF = NULL;

				if ((tempCF = insert(p, curnode->child[insertpos], insertpos, newtimestamp)) == NULL)
				{
					if (needsSplit(curnode))
						split(curnode, whichCFofParent, newtimestamp);
					return NULL;
				}
				else
				{
					curnode->cf[insertpos]->lsx -= tempCF->lsx;
					curnode->cf[insertpos]->lsy -= tempCF->lsy;
					curnode->cf[insertpos]->ssx -= tempCF->ssx;
					curnode->cf[insertpos]->ssy -= tempCF->ssy;
					curnode->cf[insertpos]->n -= tempCF->n;
					return tempCF;
				}
			}
		}

		void insert(int x, int y, int newtimestamp)
		{
			points++;
			// TODO: update CF timestamps for features along path
			this->deltaT = (this->deltaT + newtimestamp - this->lastpointtime)/2;
			this->lastpointtime = newtimestamp;
			Point *p = new Point(x, y, newtimestamp);
			if (root == NULL)
			{
				root = new Node(NULL, m, M, l, L, true);
				root->addPoint(p);
			}
			else
			{
				Node *curnode = root;
				CF *tempCF = insert(p, curnode, curnode->size, newtimestamp);
				if (tempCF)
				{
					points--;
					delete tempCF;
				}
			}
		}

		void printTree(Node *node, int nodeno)
		{
			if (node == NULL)
				return;
			cout<<"At node "<<nodeno<<"\n";
			for(int i = 0; i < node->size; i++)
			{
				cout<<"CF "<<i<<": ";
				cout<<"n = "<<node->cf[i]->n<<" ";
				cout<<"lsx = "<<node->cf[i]->lsx<<" ";
				cout<<"lsy = "<<node->cf[i]->lsy<<" ";
				cout<<"ssx = "<<node->cf[i]->ssx<<" ";
				cout<<"ssy = "<<node->cf[i]->ssy<<"\n";
			}
			cout<<"\n";
			for(int i = 0; i < node->size; i++)
				printTree(node->child[i], ++nodeno);				
		}

		void printTree()
		{
			printTree(this->root, 0);	
		}

		CF* getRootCF()
		{
			return this->root->getCF(this->lastpointtime);
		}
};

double distance(Point p1, Point cluster, int done)
{
	return sqrt(pow(p1.x - (double)cluster.x/done,2) + pow(p1.y - (double)cluster.y/done, 2));
}

double distance(Point p, CF cf)
{
	return sqrt(pow(p.x - cf.lsx/cf.n, 2) + pow(p.y - cf.lsy/cf.n, 2));
}

int main(int argc, char *argv[])
{
	srand(time(NULL));
	timeval t, tstart;
	gettimeofday(&tstart, NULL);
	
	int numTasks;
	int rank;

	MPI_Status status;
	
	MPI::Init();
	numTasks = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
	
	if (rank == MASTER)
	{
		vector<Point> initialPoints;
		for(int i = 0; i < INITIALSIZE; i++)
		{
			int tempx, tempy;
			gettimeofday(&t, NULL);
			cin>>tempx>>tempy;
			int timestamp = (t.tv_sec - tstart.tv_sec) * 1000 + (((double) (t.tv_usec - tstart.tv_usec)) / 1000);
			Point p = Point(tempx, tempy, timestamp);
			initialPoints.push_back(p);
		}

		bool assigned[INITIALSIZE];
		memset(assigned, false, sizeof(assigned));
		assigned[0] = true;
		int done = 1;

		int startpoint = 0;
		double largest = 0;
		for(int i = 0; i < INITIALSIZE; i++)
		{
			double curdist = 0;
			for(int j = 0; j < INITIALSIZE; j++)
				curdist += distance(initialPoints[i], initialPoints[j], 1);
			if (curdist > largest)
			{
				largest = curdist;
				startpoint = i;
			}
		}

		Point curcenter = initialPoints[startpoint];
		while(done < INITIALSIZE/2)
		{
			double curmin = distance(initialPoints[1], curcenter, done);
			int minpos = 1;
			for(int i = 2; i < INITIALSIZE; i++)
			{
				if (distance(initialPoints[i], curcenter, done) < curmin)
				{
					curmin = distance(initialPoints[i], curcenter, done);
					minpos = i;
				}
			}
			assigned[minpos] = true;
			done++;
			curcenter = Point(curcenter.x + initialPoints[minpos].x, curcenter.y + initialPoints[minpos].y, 0);
		}
		
		vector<Point> group1, group2;
		double lsx1, lsx2, lsy1, lsy2, ssx1, ssx2, ssy1, ssy2;
		for(int i = 0; i < INITIALSIZE; i++)
		{
			if(assigned[i])
			{
				group1.push_back(initialPoints[i]);
				lsx1 += initialPoints[i].x;
				lsy1 += initialPoints[i].y;
				ssx1 += initialPoints[i].x * initialPoints[i].x;
				ssy1 += initialPoints[i].y * initialPoints[i].y;
			}
			else
			{
				group2.push_back(initialPoints[i]);
				lsx2 += initialPoints[i].x;
				lsy2 += initialPoints[i].y;
				ssx2 += initialPoints[i].x * initialPoints[i].x;
				ssy2 += initialPoints[i].y * initialPoints[i].y;
			}
		}
		
		gettimeofday(&t, NULL);
		int timestamp = (t.tv_sec - tstart.tv_sec) * 1000 + (((double) (t.tv_usec - tstart.tv_usec)) / 1000);
		CF root1 = CF(INITIALSIZE/2, lsx1, lsy1, ssx1, ssy1, timestamp);
		CF root2 = CF(INITIALSIZE/2, lsx2, lsy2, ssx2, ssy2, timestamp);
		
		int countsent1 = INITIALSIZE/2, countsent2 = INITIALSIZE/2;
		double buffer1[20], buffer2[20];
		MPI::Request r1, r2;
		
		// Write data to corresponding buffers
		buffer1[0] = group1[0].x;
		buffer2[0] = group2[0].x;
		buffer1[1] = group1[0].y;
		buffer2[1] = group2[0].y;
		buffer1[2] = group1[0].timestamp;
		buffer2[2] = group2[0].timestamp;

		r1 = MPI::COMM_WORLD.Isend(buffer1, 3, MPI::DOUBLE, 1, 0);
		r2 = MPI::COMM_WORLD.Isend(buffer2, 3, MPI::DOUBLE, 2, 0);
		for (int i = 1, j = 1; i < group1.size() || j < group2.size(); i++, j++)
		{
			if (i < group1.size())
			{
				r1.Wait();
				buffer1[0] = group1[i].x;
				buffer1[1] = group1[i].y;
				buffer1[2] = group1[i].timestamp;
				r1 = MPI::COMM_WORLD.Isend(buffer1, 3, MPI::DOUBLE, 1, 0);
			}
			if (j < group2.size())
			{
				r2.Wait();
				buffer2[0] = group2[i].x;
				buffer2[1] = group2[i].y;
				buffer2[2] = group2[i].timestamp;
				r2 = MPI::COMM_WORLD.Isend(buffer2, 3, MPI::DOUBLE, 2, 0);
			}
		}
		
		do
		{
			int x, y;
			cin>>x>>y;
			if (cin.eof())
				break;
			gettimeofday(&t, NULL);
			int timestamp = (t.tv_sec - tstart.tv_sec) * 1000 + (((double) (t.tv_usec - tstart.tv_usec)) / 1000);
			Point p(x, y, timestamp);
			if (countsent1 == UPDATECHECK)
			{
				r1.Wait();
				MPI::COMM_WORLD.Recv(buffer1, 5, MPI::DOUBLE, 1, 0);
				root1 = CF(buffer1[0], buffer1[1], buffer1[2], buffer1[3], buffer1[4], timestamp); 
			}
			if (countsent2 == UPDATECHECK)
			{
				r2.Wait();
				MPI::COMM_WORLD.Recv(buffer2, 5, MPI::DOUBLE, 2, 0);
				root2 = CF(buffer2[0], buffer2[1], buffer2[2], buffer2[3], buffer2[4], timestamp);
			}
			if (distance(p, root1) < distance(p, root2))
			{
				r1.Wait();
				buffer1[0] = p.x;
				buffer1[1] = p.y;
				buffer1[2] = p.timestamp;
				r1 = MPI::COMM_WORLD.Isend(buffer1, 3, MPI::DOUBLE, 1, 0);
			}
			else
			{
				r2.Wait();
				buffer2[0] = p.x;
				buffer2[1] = p.y;
				buffer2[2] = p.timestamp;
				r2 = MPI::COMM_WORLD.Isend(buffer2, 3, MPI::DOUBLE, 2, 0);
			}
		} while(true);
		
		// End master node's work
	}
	else if (rank == SLAVE1)
	{
		double buffer[20];
		int countreceived = 5;
		ClusTree T(1, 3, 1, 3);
		do
		{
			MPI::COMM_WORLD.Recv(buffer, 3, MPI::DOUBLE, 0, 0);
			double x = buffer[0], y = buffer[1], timestamp = buffer[2];
			cout<<"Slave 1 received "<<x<<" "<<y<<" "<<timestamp<<"\n";
			T.insert(x, y, timestamp);
			countreceived++;
			if (countreceived == UPDATECHECK)
			{
				countreceived = 0;
				CF *tempCF = T.getRootCF();
				buffer[0] = tempCF->n;
				buffer[1] = tempCF->lsx;
				buffer[2] = tempCF->lsy;
				buffer[3] = tempCF->ssx;
				buffer[4] = tempCF->ssy;
				delete tempCF;
				MPI::COMM_WORLD.Send(buffer, 5, MPI::DOUBLE, 0, 1);
			}
		} while (true);

		// End slave 1 code
	}
	else if (rank == SLAVE2)
	{
		double buffer[20];
		int countreceived = 5;
		ClusTree T(1, 3, 1, 3);
		do
		{
			MPI::COMM_WORLD.Recv(buffer, 3, MPI::DOUBLE, 0, 0);
			double x = buffer[0], y = buffer[1], timestamp = buffer[2];
			cout<<"Slave 2 received "<<x<<" "<<y<<" "<<timestamp<<"\n";
			T.insert(x, y, timestamp);
			countreceived++;
			if (countreceived == UPDATECHECK)
			{
				countreceived = 0;
				CF *tempCF = T.getRootCF();
				buffer[0] = tempCF->n;
				buffer[1] = tempCF->lsx;
				buffer[2] = tempCF->lsy;
				buffer[3] = tempCF->ssx;
				buffer[4] = tempCF->ssy;
				delete tempCF;
				MPI::COMM_WORLD.Send(buffer, 5, MPI::DOUBLE, 0, 2);
			}
		} while (true);

		// End slave 2 code
	}

	MPI::Finalize();
	return 0;		
}
