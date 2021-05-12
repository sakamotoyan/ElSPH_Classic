#pragma once

#include "elGlobalValue.h"

/****************** DATA BLOCK ******************/
// DataBlock class stores a 2D matrix
// Each column represents a particle
// Eeach row represents one element of a specific particle
// Sereval objects of DataBlock can be assembled into one DataGrid
template <class T>
class DataBlock
{
	typedef Array<T, Dynamic, Dynamic> ArrayXXT;
	typedef Array<T, Dynamic, 1> ArrayXT;

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	/*** Section 1 all datas ***/
	ArrayXXT dataMat;

	/*** Section 2 utility functions ***/
	inline const val_i getNumElements() { return dataMat.rows(); }
	inline const val_i getNumMembers() { return dataMat.cols(); }

	inline void setConstant(T value) { dataMat.setConstant(value); }

	const inline val_i index(Array3_i& blockPseudoStruct, val_i i, val_i j, val_i k)
	{
		return i + (blockPseudoStruct(POS_X) * j) + (blockPseudoStruct(POS_X) * blockPseudoStruct(POS_Y) * k);
	}

	inline void expand(val_i rc, val_i num, val_f value = 0)
	{
		expandMat(dataMat, rc, num, value);
	}

	inline void contract(val_i rc, val_i num)
	{
		contractMat(dataMat, rc, num);
	}

	inline void show(val_i startRow = 0, val_i startCol = 0, val_i numRows = -1, val_i numCols = -1)
	{
		if (numRows == -1)
			numRows = dataMat.rows();
		if (numCols == -1)
			numCols = dataMat.cols();
		cout << dataMat.block(startRow, startCol, numRows, numCols) << endl;
	}

	/*** Section 3 constructors (and ohter means of construction) ***/

	// 3.1
	inline DataBlock() {};

	// 3.2
	inline DataBlock(val_i e, val_i a)
	{
		dataMat.resize(e, a);
	};

	// 3.3
	inline DataBlock(const DataBlock& dataBlock2)
	{
		dataMat.resize(dataBlock2.getNumElements(), dataBlock2.getNumMembers());
		dataMat = dataBlock2.dataMat;
	};

	// 3.4 positionBlock

	inline void asPosBlock(const Array3_f& cellSize, const Array3_f& centerPos, const Array3_f& cuboidSize, bool ceil = false)
	{
		val_i numElements = 3;
		Array3_i blockPseudoStruct;
		if (ceil) {
			blockPseudoStruct = (cuboidSize / cellSize).ceil().cast<val_i>();
		}
		else {
			blockPseudoStruct = (cuboidSize / cellSize).cast<val_i>();
		}
		dataMat.resize(numElements, iterProduct(blockPseudoStruct));

		Array3_f deviationPos((centerPos - (blockPseudoStruct.cast<val_f>() * cellSize / 2) + (cellSize / 2)));

#pragma omp parallel for
		for (val_i z = 0; z < blockPseudoStruct(POS_Z); z++) {
			for (val_i y = 0; y < blockPseudoStruct(POS_Y); y++) {
				for (val_i x = 0; x < blockPseudoStruct(POS_X); x++) {
					dataMat.col(index(blockPseudoStruct, x, y, z)) = deviationPos + cellSize * Array3_f(x, y, z);
				}
			}
		}
	}

	// 3.5 elementBlock
	template <class BlockT>
	inline void asElementBlock(val_i elementNum, DataBlock<BlockT>& refBlock, BlockT value = 0)
	{
		dataMat.resize(elementNum, refBlock.getNumMembers());
		dataMat.setConstant(value);
	}

	/*** Section 4 block operations ***/

	inline void resize(val_i row, val_i col)
	{
		dataMat.resize(row, col);
	}

	inline void addElement(val_i num, T value = 0)
	{
		expandMat(dataMat, ROW, num, value);
	}

	inline void pushBlock(DataBlock& dataBlock2)
	{
		val_i blockNumbersPushed = dataBlock2.getNumMembers();
		if (this->getNumElements() != dataBlock2.getNumElements()) {
			cout << "Datastruct push() error, element number didn't match" << endl;
			exit(0);
		}
		dataMat.conservativeResize(this->getNumElements(), this->getNumMembers() + blockNumbersPushed);
		dataMat.block(0, this->getNumMembers() - blockNumbersPushed, this->getNumElements(), blockNumbersPushed) =
			dataBlock2.dataMat.block(0, 0, this->getNumElements(), blockNumbersPushed);
	}

	inline void pullNumber(ArrayXi& rmNumberList)
	{
		val_i numRemoved = rmNumberList.rows();
		val_i numOrigin = dataMat.cols();
		val_i timer = 0;
		ArrayXXT tmpMat;
		Array<val_i, Dynamic, Dynamic> surviveList;
		surviveList.resize(numOrigin, 1);

		tmpMat.resize(dataMat.rows(), numOrigin - numRemoved);
		surviveList.setConstant(true);

#pragma omp parallel for
		for (val_i i = 0; i < numRemoved; i++) {
			surviveList(rmNumberList(i)) = -1;
		}

		// cannot be paralleled
		for (val_i i = 0; i < numOrigin; i++) {
			if (surviveList(i) == -1) { timer++; }
			else { surviveList(i) = i - timer; }
		}

#pragma omp parallel for
		for (val_i i = 0; i < numOrigin; i++) {
			if (surviveList(i) > -1)
				tmpMat.col(surviveList(i)) = dataMat.col(i);
		}

		dataMat = tmpMat;
	};
};

typedef DataBlock<val_f> FloatBlock;
typedef DataBlock<val_i> IntBlock;
typedef DataBlock<double> DataBlock_d;
typedef DataBlock<float> DataBlock_f;
typedef DataBlock<val_i> DataBlock_i;

/****************** DATA GRID ******************/
class DataGrid
{
public:
	/*** Section 1 all datas ***/

	vector<IntBlock*> intBlocks;
	vector<FloatBlock*> floatBlocks;

	val_i configInfo[3];

	// config of the marker blocks
	IntBlock infoBlock;		// 1st marker block
	// config of the value blocks
	FloatBlock SoABlock;			// 1st value block
	FloatBlock AoSBlock;			// 2nd value block

	// upper limit
	val_i numPart = 0;

	/*** Section 2 constructors (and ohter means of construction) ***/
	inline DataGrid(val_i num = 0, const val_i* config = basicGridConfig)
	{
		numPart = num;
		configInfo[0] = config[0];
		configInfo[1] = config[1];
		configInfo[2] = config[2];

		intBlocks.push_back(&infoBlock);
		floatBlocks.push_back(&SoABlock);
		floatBlocks.push_back(&AoSBlock);

		for (val_i i = 0; i < intBlocks.size(); i++) {
			intBlocks[i]->resize(config[i], num);
			intBlocks[i]->dataMat.setZero();
		}

		for (val_i i = 0; i < floatBlocks.size(); i++) {
			floatBlocks[i]->resize(config[i + intBlocks.size()], num);
			floatBlocks[i]->dataMat.setZero();
		}
	};

	inline void reset(val_i num = 0) {
		for (val_i i = 0; i < intBlocks.size(); i++) {
			intBlocks[i]->resize(configInfo[i], num);
			intBlocks[i]->dataMat.setZero();
		}

		for (val_i i = 0; i < floatBlocks.size(); i++) {
			floatBlocks[i]->resize(configInfo[i + intBlocks.size()], num);
			floatBlocks[i]->dataMat.setZero();
		}
	}

	inline void clean() {
		for (val_i i = 0; i < intBlocks.size(); i++) {
			intBlocks[i]->dataMat.setZero();
		}

		for (val_i i = 0; i < floatBlocks.size(); i++) {
			floatBlocks[i]->dataMat.setZero();
		}
	}

	inline void clean(val_i colNum) {
		for (val_i i = 0; i < floatBlocks.size(); i++) {
			floatBlocks[i]->dataMat.block(0, 0, floatBlocks[i]->dataMat.rows(), colNum).setZero();
		}
		for (val_i i = 0; i < intBlocks.size(); i++) {
			intBlocks[i]->dataMat.block(0, 0, intBlocks[i]->dataMat.rows(), colNum).setZero();
		}
	}

	inline void cleanSpecificFloatAttributes(val_i* attrRef) {
		floatBlocks[attrRef[0]]->dataMat.block(attrRef[1], 0, attrRef[2], numPart).setZero();
	}
	inline void cleanSpecificIntAttributes(val_i* attrRef) {
		intBlocks[attrRef[0]]->dataMat.block(attrRef[1], 0, attrRef[2], numPart).setZero();
	}

	inline const val_i num() { return SoABlock.getNumMembers(); }

	/*** Section 4 grid operations ***/
	inline void pushGrid(DataGrid& pushedGrid) {
		for (val_i i = 0; i < intBlocks.size(); i++) {
			intBlocks[i]->pushBlock(*(pushedGrid.intBlocks[i]));
		}
		for (val_i i = 0; i < floatBlocks.size(); i++) {
			floatBlocks[i]->pushBlock(*(pushedGrid.floatBlocks[i]));
		}
	}

	inline void setConstantVal(val_f value) {
		for (val_i i = 0; i < floatBlocks.size(); i++) {
			floatBlocks[i]->setConstant(value);
		}
	}

	inline void setConstantMarker(val_i value) {
		for (val_i i = 0; i < intBlocks.size(); i++) {
			intBlocks[i]->setConstant(value);
		}
	}

	/*** Section 3 utility functions ***/

	inline void show(val_i showCase = -1)
	{
		cout << "/-----------------------------" << endl;
		for (val_i i = 0; i < intBlocks.size(); i++) {
			cout << "|*** integer block [" << i << "] ***" << endl;
			intBlocks[i]->show(0, 0, -1, showCase);
		}
		cout << "-----------------------------" << endl;
		for (val_i i = 0; i < floatBlocks.size(); i++) {
			cout << "*** floating block [" << i << "] ***" << endl;
			floatBlocks[i]->show(0, 0, -1, showCase);
		}
		cout << "|" << endl << "\\-----------------------------" << endl;
		cout << endl;
	}
};


