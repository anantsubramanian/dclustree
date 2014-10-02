#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>

using namespace std;

#define LAMBDA 0.00001
#define BETA 2
#define NUMINSERTS 300

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

			// TODO: Check timestamp
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
};

int main()
{
	srand(time(NULL));
	timeval t, tstart;
	ClusTree T(1, 3, 1, 3);
	int x, y, timestamp, count = 0;
	gettimeofday(&tstart, NULL);
	int starttime, endtime;
	bool first = true;
	do
	{
		cin>>x>>y;
		if (cin.eof())
			break;
		gettimeofday(&t, NULL);
		int timestamp = (t.tv_sec - tstart.tv_sec) * 1000 + (((double) (t.tv_usec - tstart.tv_usec)) / 1000);
		if (first)
		{
			first = false;
			starttime = timestamp;
		}
		T.insert(x, y, timestamp);
		points++;
		//cout<<points<<"\n";
	} while(true);
	cout<<points<<"\n";
	//gettimeofday(&t, NULL);
	//timestamp = (t.tv_sec - tstart.tv_sec) * 1000 + (((double) (t.tv_usec - tstart.tv_usec)) / 1000);
	//endtime = timestamp;
	//cout<<endtime - starttime<<"\n";
	return 0;		
}
