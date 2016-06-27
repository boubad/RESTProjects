#include "stdafx.h"
#include "CppUnitTest.h"
////////////////////////////////////
#include "../geninfo/include/matdata.h"
//////////////////////////////////
#include "infotestdata.h"
////////////////////////////////////
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace info;
namespace UnitTestGenInfo
{		
	//
	using string_type = utility::string_t;
	using strings_vector = std::vector<string_type>;
	using doubles_vector = std::vector<double>;
	using DistanceMapType = DistanceMap<double>;
	using PDistanceMapType = DistanceMapType *;
	using MatDataPtr = std::shared_ptr<MatData>;
	//
	TEST_CLASS(UnitTestMatData)
	{
	public:
		
		TEST_METHOD(MatDataInitAsync)
		{
			string_type name;
			size_t nRows = 0, nCols = 0;
			std::vector<int> data;
			strings_vector rowNames, colNames;
			//
			InfoTestData::get_mortal_data(name, nRows, nCols, data, rowNames, colNames);
			//
			pplx::task<MatDataPtr> tsk = MatData::createAsync(nRows, nCols, data);
			MatDataPtr oMat = tsk.get();
			MatData *pData = oMat.get();
			Assert::IsNotNull(pData);
			pData->set_names(&rowNames, &colNames);
		}// MatDataInitAsync

	};
}