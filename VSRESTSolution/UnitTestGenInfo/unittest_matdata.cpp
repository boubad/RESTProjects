#include "stdafx.h"
#include "CppUnitTest.h"
////////////////////////////////////
#include "../geninfo/include/matdata.h"
#include "../geninfo/include/matelem.h"
//////////////////////////////////
#include "infotestdata.h"
#include "global_defs.h"
////////////////////////////////////
using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace info;
///////////////////////////////////////////////
namespace UnitTestGenInfo
{	
	//
	using MatDataType = MatData<INDEXTYPE,FLOATTYPE, DISTANCETYPE>;
	//
	using string_type = typename MatDataType::string_type;
	using strings_vector = typename MatDataType::strings_vector;
	using DistanceMapType = typename MatDataType::DistanceMapType;
	using PDistanceMapType = typename MatDataType::PDistanceMapType;
	using MatDataPtr = typename MatDataType::MatDataPtr;
	using MatArgsType = MatArgs<DATATYPE>;
	//
	using MatElemType = MatElem<DISTANCETYPE, CRITERIATYPE>;
	using index_ptr_type = typename MatElemType::index_ptr_type;
	using result_type = typename MatElemType::result_type;
	using task_type = pplx::task<result_type>;
	using criteria_type = typename MatElemType::criteria_type;
	using sizets_vector = typename MatElemType::sizets_vector;
	//
	TEST_CLASS(UnitTestMatData)
	{
	public:
		
		TEST_METHOD(MatDataInitAsync)
		{
			string_type name;
			size_t nRows = 0, nCols = 0;
			std::vector<DATATYPE> data;
			strings_vector rowNames, colNames;
			//
			InfoTestData::get_mortal_data(name, nRows, nCols, data, rowNames, colNames);
			std::shared_ptr<MatArgsType> oArgs = std::make_shared<MatArgsType>(nRows, nCols, data, &rowNames, &colNames);
			//
			MatDataPtr oMat  = MatDataType::create(oArgs);
			MatDataType *pData = oMat.get();
			Assert::IsNotNull(pData);
			PDistanceMapType pMap = pData->get_rows_distances_map();
			Assert::IsNotNull(pMap);
			MatElemType xMat(pMap);
			//
			result_type oRes = xMat.arrange();
			criteria_type c = oRes.first;
			Assert::IsTrue(c > 0);
			index_ptr_type oind = oRes.second;
			const sizets_vector *ps = oind.get();
			Assert::IsNotNull(ps);
		}// MatDataInitAsync

	};
}