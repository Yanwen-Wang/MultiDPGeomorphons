/***************************************************************************
* pAspect.cpp
*
* Project: GPRO, v 1.0
* Purpose: Demonstration program for GPRO. 
*
* Usage:  
*
* Example: 
*
* Author:  Wang Yanwen
* E-mail:  wangyw@lreis.ac.cn
****************************************************************************
* Copyright (c) 2015. Zhan Lijun
* 
****************************************************************************/


#include <iostream>
#include <string>
#include <fstream>
#include <time.h>
#include <sstream>
#include<omp.h>
#include "mpi.h"
#include "neighborhood.h"
#include "cellSpace.h"
#include "basicTypes.h"
#include "rasterLayer.h"
#include "application.h"
#include "multiDPGMPOperator.h"
#include "communication.h"


using namespace std;
using namespace GPRO;

int main(int argc, char *argv[]) 
{
	/*  enum ProgramType{MPI_Type = 0,
				   MPI_OpenMP_Type,
				   CUDA_Type,
				   Serial_Type};*/
	Application::START(MPI_Type, argc, argv); //init

	//...
	char* inputfilename;
	char* neighborfile;
	char* outputfilename;
	//char* outputfile2name;
	int threadNUM;
	if (argc < 6)
	{
		inputfilename = argv[1];
		neighborfile = argv[2]; 
		outputfilename = argv[3];
		//outputfile2name = argv[4];
	}
	//omp_set_num_threads(threadNUM);
	RasterLayer<double> demLayer("demLayer"); //创建图层
	demLayer.readNeighborhood(neighborfile);  //读取分析窗口文件
	demLayer.readFile(inputfilename);  //读取栅格数据

	RasterLayer<double> slopeLayer("slopeLayer");
	//RasterLayer<double> AcurracyLayer("AcurracyLayer");

	slopeLayer.copyLayerInfo(demLayer);
	//AcurracyLayer.copyLayerInfo(demLayer);

	double starttime;
	double endtime;

	MPI_Barrier(MPI_COMM_WORLD);
	starttime = MPI_Wtime();
	multiDPGMPOperator SlpOper;
	
	//cout<<"SlopeOperator SlpOper;"<<endl;
	SlpOper.demLayer(demLayer);
	
	SlpOper.slopeLayer(slopeLayer);
	
	//cout<<"SlpOper.slopeLayer(slopeLayer);"<<endl;
	SlpOper.Run();
	
	MPI_Barrier(MPI_COMM_WORLD);

	endtime = MPI_Wtime();


	cout<<"run time is "<<endtime-starttime<<endl;
	slopeLayer.writeFile(outputfilename);
	
	
	Application::END();
	return 0;
}
