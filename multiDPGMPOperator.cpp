  #include"multiDPGMPOperator.h"
#include"math.h"
#include"geomorphons.h"
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>

//为了控制DP的递归次数
/*directions
 * 3|2|1
 * 4|0|8
 * 5|6|7 */
//int nextr[8] = {-1, -1, -1, 0, 1, 1, 1, 0 };
//int nextc[8] = {1, 0, -1, -1, -1, 0, 1, 1 };

void multiDPGMPOperator::demLayer(RasterLayer<double> &layerD) 
{
  _pDEMLayer = &layerD;
  _pDEMNbrhood = layerD.nbrhood();
  cellSize = _pDEMLayer->_pMetaData->cellSize;
  noData = _pDEMLayer->_pMetaData->noData;
    //pWorkBR = &_pDEMLayer->_pMetaData->_localworkBR;
  Configure(_pDEMLayer, false);
}

void multiDPGMPOperator::slopeLayer(RasterLayer<double> &layerD) 
{
  _pSlopeLayer = &layerD;
  Configure(_pSlopeLayer,false);
}

bool multiDPGMPOperator::isTermination()
{
	num--;
	if(num > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

//计算距离
double caldist(int row1,int col1,int row2,int col2,double cell_res)
{
	double dist;
	dist = cell_res * sqrt((double)((col1-col2)*(col1-col2) +(row1-row2)*(row1-row2)));
	return dist;
}
//计算线参数
void ParaCalcu(vector<GeoPoint> Points,int firstpoint, int lastpoint,double *k,double *b)
{
	//LineParameters thisLine = new LineParameters();
	*k = ((double)Points[lastpoint].elevation-(double)Points[firstpoint].elevation)/((double)Points[lastpoint].distance-(double)Points[firstpoint].distance);
	*b = (double)Points[firstpoint].elevation - *k*(double)Points[firstpoint].distance;
	//cout<<"step1.1.1.1"<<" "<<firstpoint<<" "<<lastpoint<<endl;
	//cout<<"step1.1.1"<<" "<<Points[lastpoint].elevation<<" "<<Points[lastpoint].distance<<" "<<Points[firstpoint].elevation<<" "<<Points[firstpoint].distance<<endl;
	//cout<<"step1.1"<<" "<<k<<" "<<b<<endl;
}
//计算每一点离直线的高差
double DisCalcu(vector<GeoPoint> Points, int thispoint, double k,double b)
{
	//GeoPoint temPoint = new GeoPoint(0,0,0);
	double tem = 0; double distance;			
	tem = Points[thispoint].distance * k + b;			
	//distance = fabs(Points[thispoint].elevation - tem);
	double up = fabs((-k)*Points[thispoint].distance + Points[thispoint].elevation - b);
	//cout<<"up="<<up<<endl;
	double down = sqrt(k*k+1);
	//cout<<"down="<<down<<endl;
	distance = up/down;
	//cout<<"distance="<<distance<<endl;
	return distance;
}
//进行DP一次计算
void DPCalcu(vector<GeoPoint> Points, int firstpoint, int lastpoint, double Tolerance, vector<int> &PointToKeep)
{
	double maxDistance = 0;
	int FarthestPoint = 0;
	double k,b;

	//计算每一点离直线的高差
	ParaCalcu(Points, firstpoint, lastpoint,&k,&b);
	//cout<<"step1"<<k<<" "<<b<<endl;
	for (int thispoint = firstpoint;thispoint<lastpoint;thispoint++)
	{
		double distance = DisCalcu(Points,thispoint,k,b);
		//cout<<"step2"<<thispoint<<" "<<distance<<endl;
		if (distance > maxDistance)
		{
			maxDistance = distance;
			FarthestPoint = thispoint;
		}
		//cout<<"step3"<<maxDistance<<" "<<FarthestPoint<<endl;
	}

		//将最大距离点保留
		PointToKeep.push_back(FarthestPoint);
		//递归计算
		//if (times == 0)
		//{
			//times = times + 1;
			//DPCalcu(Points, firstpoint, FarthestPoint, Tolerance,PointToKeep,times);
			//DPCalcu(Points, FarthestPoint, lastpoint, Tolerance, PointToKeep,times);
		//}
	//else
		return;
}

//化简patten类型
int determine_form(int num_minus, int num_plus)
{
/* determine form according number of positives and negatives
 * simple approach to determine form pattern */
	const FORMS forms[9][9] = 
	{
/* minus ------------- plus ----------------*/
/*       0   1   2   3   4   5   6   7   8  */
/* 0 */ {FL, FL, FL, FS, FS, VL, VL, VL, PT},
/* 1 */ {FL, FL, FS, FS, FS, VL, VL, VL, __},
/* 2 */ {FL, SH, SL, SL, CN, CN, VL, __, __},
/* 3 */ {SH, SH, SL, SL, SL, CN, __, __, __},
/* 4 */ {SH, SH, CV, SL, SL, __, __, __, __},
/* 5 */ {RI, RI, CV, CV, __, __, __, __, __},
/* 6 */ {RI, RI, RI, __, __, __, __, __, __},
/* 7 */ {RI, RI, __, __, __, __, __, __, __},
/* 8 */ {PK, __, __, __, __, __, __, __, __},
    };

/* legend:
  FL,  flat
  PK,  peak, summit
  RI,  ridge
  SH,  shoulder
  CV,  convex
  SL,  slope
  CN,  concave
  FS,  footslope
  VL,  valley
  PT,  pit, depression
  __  error, impossible
*/
 return forms[num_minus][num_plus];
}

int determine_form_num(int num_minus, int num_plus)
{
/* determine form according number of positives and negatives
 * simple approach to determine form pattern */
	const FORMS_NUM forms_num[9][9] = 
	{
		/* minus ------------- plus ----------------*/
		/*       0   1   2   3   4   5   6   7   8  */
		/* 0 */ {FLN, FLN, FLN, FSN, FSN, VLN, VLN, VLN, PTN},
		/* 1 */ {FLN, FLN, FSN, FSN, FSN, VLN, VLN, VLN, __N},
		/* 2 */ {FLN, SHN, SLN, SLN, CNN, CNN, VLN, __N, __N},
		/* 3 */ {SHN, SHN, SLN, SLN, SLN, CNN, __N, __N, __N},
		/* 4 */ {SHN, SHN, CVN, SLN, SLN, __N, __N, __N, __N},
		/* 5 */ {RIN, RIN, CVN, CVN, __N, __N, __N, __N, __N},
		/* 6 */ {RIN, RIN, RIN, __N, __N, __N, __N, __N, __N},
		/* 7 */ {RIN, RIN, __N, __N, __N, __N, __N, __N, __N},
		/* 8 */ {PKN, __N, __N, __N, __N, __N, __N, __N, __N},
	};

/* legend:
  FL,  flat
  PK,  peak, summit
  RI,  ridge
  SH,  shoulder
  CV,  convex
  SL,  slope
  CN,  concave
  FS,  footslope
  VL,  valley
  PT,  pit, depression
  __  error, impossible
*/
 return forms_num[num_minus][num_plus];
}

//三窗口计算form
int extern multidetermin_form_num(Pattern patterns[],int RowNum)
{
	int forms[RowNum-2];
	//int finalform = 0;
	
	for (int i=0;i<RowNum-2;i++)
	{
		forms[i] = determine_form_num(patterns[i].num_negatives,patterns[i].num_positives);
	}
	//当三个窗口都一样或者都不一样时，选择forms[0]为finalform
	//if (forms[0]==forms[1]==forms[2])
	//{
		//finalform = forms[0];
	//}
	//if (forms[0]!=forms[1]&&forms[0]!=forms[2]&&forms[1]!=forms[2])
	//{
		//finalform = forms[0];
	//}

	//记录各种form出现次数,算上错误类型一共11种类型
	int curform[11] = {0,0,0,0,0,0,0,0,0,0,0};
	//对三个窗口的form进行计算
	
	for (int i = 0;i<RowNum-2;i++)
	{
		switch(forms[i])
		{
		case 1:curform[0]++;break;
		case 2:curform[1]++;break;
		case 3:curform[2]++;break;
		case 4:curform[3]++;break;
		case 5:curform[4]++;break;
		case 6:curform[5]++;break;
		case 7:curform[6]++;break;
		case 8:curform[7]++;break;
		case 9:curform[8]++;break;
		case 10:curform[9]++;break;
		default:curform[10]++;break;
		}
	}

	int formmax = 0;
	int formmaxnum = 0;
	for (int i=0;i<11;i++)
	{
		if (curform[i]>formmax)
		{
			formmax = i+1;
			formmaxnum = curform[i];
		}
	}
	int form = formmax;
	//int form = formmaxnum;
	//cout<<"form = "<<form<<endl;
	//form[1] = curform[i];
	return form;
}

bool multiDPGMPOperator::Operator(const CellCoord &coord)
{
	CellSpace<double> &dem = *(_pDEMLayer->cellSpace());//输入图层的栅格数据
	
	CellSpace<double> &slope = *(_pSlopeLayer->cellSpace());//输出图层的栅格数据
	  
	Neighborhood<double>& nbrhoodD = *(_pDEMNbrhood);//分析窗口文件

	int nextr[8] = {-1, -1, -1, 0, 0, 1, 1, 1 };
	int nextc[8] = {-1, 0, 1, -1, 1, -1, 0, 1 };
	
	//coord是分析窗口的坐标原点，也就是需要计算坡度的点
	int cur_Row = coord.iRow();//当前的row
	int cur_Col = coord.iCol();//当前的cols

	//neighborhood文件一共有十个点，前九个都是步长点，第十个是窗口大小点
	int RowNum = sqrt(nbrhoodD.size()-1)/2;//窗口的row数
	int ColNum = sqrt(nbrhoodD.size()-1)/2;//窗口的col数

	//double cell_res = _pDEMLayer->_pMetaData->cellSize;//	DEM格网长度
	double cell_res = 30;
	/*
	Pattern pattern;//pattern初始化
	pattern.num_negatives = 0;
	pattern.num_positives = 0;
	pattern.negatives = NULL;
	pattern.positives = NULL;
	//初始化三窗口pattern
	Pattern pattern2,pattern3;
	pattern2 = pattern3 = pattern;
	*/

	Pattern patterns[RowNum-2];
	for (int i=0;i<RowNum-2;i++)
	{
		patterns[i].num_negatives = 0;
		patterns[i].num_positives = 0;
		patterns[i].negatives = 0;
		patterns[i].positives = 0;
	}
	//设置角度阈值
	//double flat_distance = cell_res*3;//为了避免窗口内起伏的误差，所有点将从3格之外开始计算
	//double flat_threshold = 0.017;//设置为1°
	//double flat_threshold_height = flat_distance*tan(flat_threshold);//随便设置一个数
	//将角度阈值设置为3°
	double flat_distance = cell_res*3;//为了避免窗口内起伏的误差，所有点将从3格之外开始计算
	double flat_threshold = 0.052;//设置为3°
	double flat_threshold_height = flat_distance*tan(flat_threshold);//随便设置一个数
	

	vector<GeoPoint> Points;//创造点集
	vector<int> PointToKeep;//创造保留点集

	//控制递归次数
	//int times;

	for(int i=0;i<8;i++)
	{
		//多窗口计算，因此需要从第3个点计算到第RowNum个点
		for (int j = 0;j<RowNum-2;j++)
		{
			patterns[j].pattern[i]=0;

			//初始化Points
			GeoPoint curPoint;
			curPoint.i = i;
			//每次循环开始让times等于0控制递归次数
			//times = 0;

			for (int k = 0;k<j+3;k++)
			{
				//cout<<"cell_res = "<<cell_res<<endl;
				curPoint.row = cur_Row + k*nextr[i];
				//cout<<"step2.1.1"<<" "<<k<<" "<<curPoint.row;
				curPoint.col = cur_Col + k*nextc[i];
				//cout<<"step2.1.2"<<" "<<k<<" "<<curPoint.col;
				curPoint.elevation = dem[curPoint.row][curPoint.col];
				//cout<<"step2.1.3"<<" "<<k<<" "<<curPoint.elevation<<" ";
				curPoint.height = curPoint.elevation - dem[cur_Row][cur_Col];
				//cout<<"step2.1.4"<<" "<<k<<" "<<curPoint.height;
				//curPoint.distance = caldist(curPoint.row,curPoint.col,cur_Row,cur_Col,cell_res);
				curPoint.distance = cell_res * sqrt((double)((curPoint.col-cur_Col)*(curPoint.col-cur_Col) +(curPoint.row-cur_Row)*(curPoint.row-cur_Row)));
				//cout<<"step2.1"<<" "<<k<<" "<<curPoint.distance<<endl;
				curPoint.angel = atan2(curPoint.height,curPoint.distance);
				Points.push_back(curPoint);
			}

			//初始化计算所需点
			double Tolerance = 0;
			int firstpoint = 0;
			int lastpoint  = Points.size()-1;

			//PointToKeep.push_back(firstpoint);
			//PointToKeep.push_back(lastpoint);

			DPCalcu(Points,firstpoint,lastpoint,Tolerance,PointToKeep);

			int feature_point = PointToKeep[0];
			//int feature_point2,feature_point3;
			//feature_point2 = PointToKeep[1];
			//feature_point3 = PointToKeep[2];
			//cout<<feature_point<<"  "<<Points[feature_point].angel<<"  "<<flat_threshold<<endl;
			if (Points[feature_point].angel>flat_threshold)
			{
				patterns[j].positives += i;
				patterns[j].pattern[i] = 1;
				patterns[j].num_positives++;
				patterns[j].feature[i] = Points[feature_point];
			}
			if (Points[feature_point].angel<-flat_threshold)
			{
				patterns[j].negatives += i;
				patterns[j].pattern[i] = -1;
				patterns[j].num_negatives++;
				patterns[j].feature[i] = Points[feature_point];
			}
			else
			{
				patterns[j].pattern[i] = 0;
				patterns[j].feature[i] = curPoint;
			}
			//清除原来的数据
			Points.clear();
			PointToKeep.clear();
		}
	}
	vector<GeoPoint>().swap(Points);
	vector<int>().swap(PointToKeep);
	int finalForm;
	finalForm = multidetermin_form_num(patterns,RowNum);
	//cout<<"final = "<<finalForm<<endl;
	slope[cur_Row][cur_Col] = finalForm;
	//slope[cur_Row][cur_Col] = finalForm[1];
	//delete[] finalForm;
	return true;
}
