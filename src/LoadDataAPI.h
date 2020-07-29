#ifndef _LOAD_DATA_API_H_
#define _LOAD_DATA_API_H_

#include <memory>

// Required types:
//
// enum ColumnDataType
// dataPtr
// struct ColumnData
// ColMapType

class LoadDataProviderBase2 {
protected:
	const char *name;
	const char *dataName;
	int destRows;  // ==sum(rowFilter)
	std::vector<ColumnData> *rawCols;
	ColMapType *rawColMap;
	std::vector< int > columns;
	std::vector< ColumnDataType > colTypes;
	std::vector<dataPtr> origData;  // who deallocates? TODO
	bool checkpoint;
	std::vector< std::string > *checkpointValues;

	int verbose;
	int curRecord;
	int loadCounter;
	int rowNames, colNames;
	int skipRows, skipCols;
	std::vector<std::string> naStrings;
	int srcRows;  // destRows or length(rowFilter)
	int *rowFilter;

	std::string filePath;
	std::string fileName;

	int stripeSize;
	int stripeStart;  // 0 is the first column
	int stripeEnd;
	std::vector<dataPtr> stripeData; // stripeSize * columns.size()

	bool isNA(const std::string& str)
	{
		for (auto &na1 : naStrings) {
			if (na1 == str) return true;
		}
		return false;
	}

	void requireFile(SEXP rObj1)
	{
		RObject rObj(rObj1);
		CharacterVector Rpath(rObj.slot("path"));
		if (Rpath.size() != 1)
			mxThrow("%s: you must specify exactly one file from which to read data", name);

		filePath = Rpath[0];
		auto slashPos = filePath.find_last_of("/\\");
		if (slashPos == std::string::npos) {
			fileName = filePath;
		} else {
			fileName = filePath.substr(slashPos+1);
		}
	}

	virtual void loadRowImpl(int index)=0;

public:
	const std::vector< int > &getColumns() { return columns; }
	void commonInit(SEXP rObj, const char *_name,
                  const char *_dataName, int rows,
                  std::vector<ColumnData> &_rawCols,
                  ColMapType &_rawColMap,
                  std::vector< std::string > &_checkpointValues,
                  bool useOriginalData);
	virtual int getNumVariants() { return 0; }
	bool wantCheckpoint() const { return checkpoint; }
	int getLoadCounter() const { return loadCounter; }
	bool skipRow(int rx) const { return !rowFilter? false : rowFilter[rx]; }
	virtual const char *getName()=0;
	virtual void init(SEXP rObj)=0;
	virtual void addCheckpointColumns(std::vector< std::string > &cp) {};
	void loadRow(int index)
	{
		if (!stripeData.size()) {
			stripeData.reserve(stripeSize * columns.size());
			for (int sx=0; sx < stripeSize; ++sx) {
				for (int cx=0; cx < int(columns.size()); ++cx) {
					if (colTypes[cx] == COLUMNDATA_NUMERIC) {
						stripeData.emplace_back(new double[destRows]);
					} else {
						stripeData.emplace_back(new int[destRows]);
					}
				}
			}
		}

		loadRowImpl(index);
	}
	void loadOrigRow() {
		auto rc = *rawCols;
		for (int cx=0; cx < int(columns.size()); ++cx) {
			rc[ columns[cx] ].setBorrow(origData[cx]);
		}
	}
	virtual std::unique_ptr<LoadDataProviderBase2> clone()=0;
	virtual ~LoadDataProviderBase2() {
		int stripes = stripeData.size() / columns.size();
		for (int sx=0; sx < stripes; ++sx) {
			for (int cx=0; cx < int(columns.size()); ++cx) {
				int dx = sx * columns.size() + cx;
				if (colTypes[cx] == COLUMNDATA_NUMERIC) {
					delete [] stripeData[dx].realData;
				} else {
					delete [] stripeData[dx].intData;
				}
			}
		}
		stripeData.clear();
	}
};

template <typename Derived>
class LoadDataProvider : public LoadDataProviderBase2 {
public:
	virtual std::unique_ptr<LoadDataProviderBase2> clone() {
		return std::unique_ptr<LoadDataProviderBase2>(new Derived());
	}
};

//#define OPENMX_LOAD_DATA_API_VERSION 0.17789282277226448059
//#define OPENMX_LOAD_DATA_API_VERSION 0.3091921037994325
#define OPENMX_LOAD_DATA_API_VERSION 0.5240939254872501 // this is a random number

typedef void (*AddLoadDataProviderType)(double version, int ldpbSz, LoadDataProviderBase2 *ldp);

#endif
