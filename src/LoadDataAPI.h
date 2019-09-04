#ifndef _LOAD_DATA_API_H_
#define _LOAD_DATA_API_H_

// Required types:
//
// enum ColumnDataType
// dataPtr
// struct ColumnData
// ColMapType

class LoadDataProviderBase {
protected:
	const char *name;
	const char *dataName;
	int rows;
	std::vector<ColumnData> *rawCols;
	ColMapType *rawColMap;
	std::vector< int > columns;
	std::vector< ColumnDataType > colTypes;
	std::vector<dataPtr> origData;
	bool checkpoint;
	std::vector< std::string > *checkpointValues;

	int verbose;
	int curRecord;
	int loadCounter;
	int rowNames, colNames;
	int skipRows, skipCols;
	std::vector<std::string> naStrings;

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

	void requireFile(SEXP rObj)
	{
		ProtectedSEXP Rpath(R_do_slot(rObj, Rf_install("path")));
		if (Rf_length(Rpath) != 1)
			mxThrow("%s: you must specify exactly one file from which to read data", name);

		filePath = R_CHAR(STRING_ELT(Rpath, 0));
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
			std::vector< std::string > &_checkpointValues);
	virtual int getNumVariants() { return 0; }
	bool wantCheckpoint() const { return checkpoint; }
	int getLoadCounter() const { return loadCounter; }
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
						stripeData.emplace_back(new double[rows]);
					} else {
						stripeData.emplace_back(new int[rows]);
					}
				}
			}
		}

		loadRowImpl(index);
	}
	void loadOrigRow() {
		auto rc = *rawCols;
		for (int cx=0; cx < int(columns.size()); ++cx) {
			rc[ columns[cx] ].ptr = origData[cx];
		}
	}
	virtual std::unique_ptr<LoadDataProviderBase> clone()=0;
	virtual ~LoadDataProviderBase() {
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
class LoadDataProvider : public LoadDataProviderBase {
public:
	virtual std::unique_ptr<LoadDataProviderBase> clone() {
		return std::unique_ptr<LoadDataProviderBase>(new Derived());
	}
};

#define OPENMX_LOAD_DATA_API_VERSION 0.17789282277226448059 // this is a random number

typedef void (*AddLoadDataProviderType)(double version, LoadDataProviderBase *ldp);

#endif
