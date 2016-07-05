#include "stdafx.h"
#include "CppUnitTest.h"
////////////////////////////////////
#include "../geninfo/include/matord.h"
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
	using index_type = INDEXTYPE;
	using criteria_type = CRITERIATYPE;
	using distance_type = DISTANCETYPE;
	using cancellableflag_type = CancellableObject::cancellableflag_type;
	using MatDataType = MatData<INDEXTYPE,FLOATTYPE, DISTANCETYPE,STRINGTYPE>;
	//
	using string_type = typename MatDataType::string_type;
	using strings_vector = typename MatDataType::strings_vector;
	using DistanceMapType = typename MatDataType::distancemap_type;
	using MatDataPtr = typename MatDataType::matdata_type_ptr;
	using matelemresult_type = MatElemResult<index_type, criteria_type>;
	using matelemresult_type_ptr = std::shared_ptr<matelemresult_type>;
	using matelem_type = MatElem<index_type, distance_type, criteria_type>;
	using matord_type = MatOrd<index_type, distance_type, criteria_type>;
	using matordresult_type = std::pair<matelemresult_type_ptr, matelemresult_type_ptr>;
	//
	
	//
	TEST_CLASS(UnitTestMatData)
	{
	public:
		void write_matelemresult(matelemresult_type_ptr r) {
			std::stringstream os;
			matelemresult_type *p = r.get();
			if (p != nullptr) {
				os << *p << std::endl;
			}
			std::string s = os.str();
			Logger::WriteMessage(s.c_str());
		 }//write_matelemresult
		void write_matordresult(matordresult_type r) {
			write_matelemresult(r.first);
			write_matelemresult(r.second);
		}//write_matelemresult
		TEST_METHOD(MatDataInitAsync)
		{
			string_type name;
			size_t nRows = 0, nCols = 0;
			std::vector<DATATYPE> data;
			strings_vector rowNames, colNames;
			//
			InfoTestData::get_mortal_data(name, nRows, nCols, data, rowNames, colNames);
			MatDataPtr oData = MatDataType::create(nRows, nCols, data, &rowNames, &colNames);
			MatDataType *pData = oData.get();
			Assert::IsNotNull(pData);
			Assert::IsTrue(pData->is_valid());
			//
		}// MatDataInitAsync
		TEST_METHOD(MatOrd)
		{
			string_type name;
			size_t nRows = 0, nCols = 0;
			std::vector<DATATYPE> data;
			strings_vector rowNames, colNames;
			//
			InfoTestData::get_mortal_data(name, nRows, nCols, data, rowNames, colNames);
			MatDataPtr oData = MatDataType::create(nRows, nCols, data, &rowNames, &colNames);
			MatDataType *pData = oData.get();
			Assert::IsNotNull(pData);
			Assert::IsTrue(pData->is_valid());
			//
			cancellableflag_type cancel(false);
			Backgrounder back;
			//
			matord_type oMat(pData,&cancel,&back);
			matordresult_type r = oMat.arrange_all();
			write_matordresult(r);
			//
		}// MatOrd
		TEST_METHOD(MatOrdHierar)
		{
			string_type name;
			size_t nRows = 0, nCols = 0;
			std::vector<DATATYPE> data;
			strings_vector rowNames, colNames;
			//
			InfoTestData::get_mortal_data(name, nRows, nCols, data, rowNames, colNames);
			MatDataPtr oData = MatDataType::create(nRows, nCols, data, &rowNames, &colNames);
			MatDataType *pData = oData.get();
			Assert::IsNotNull(pData);
			Assert::IsTrue(pData->is_valid());
			//
			matordresult_type r = matord_type::st_arrange_all_hierar(pData);
			write_matordresult(r);
			//
		}// MatOrd

	};
}