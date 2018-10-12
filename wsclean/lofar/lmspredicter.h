#ifndef LOFAR_MS_PREDICTER
#define LOFAR_MS_PREDICTER

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include "../banddata.h"
#include "../dftpredictionalgorithm.h"
#include "../lane.h"
#include "../buffered_lane.h"

#include "lbeamevaluator.h"

#include <boost/asio/io_service.hpp>
#include <boost/thread/barrier.hpp>

#include <complex>
#include <memory>
#include <thread>

class LMSPredicter
{
public:
	struct RowData
	{
		RowData() : modelData(0)
		{	}
		
		MC2x2* modelData;
		size_t rowIndex, timeIndex, a1, a2;
		double u, v, w;
	};
	
	explicit LMSPredicter(casacore::MeasurementSet &ms, size_t threadCount) :
		_ms(ms),
		_applyBeam(false),
		_useModelColumn(false),
		_dftInput(),
		_barrier(threadCount+1),
		_laneSize(threadCount*16),
		_workLane(_laneSize),
		_outputLane(_laneSize),
		_availableBufferLane(_laneSize),
		_bufferedOutputLane(&_outputLane, threadCount*2),
		_bufferedBufferLane(&_availableBufferLane, threadCount*2),
		_startRow(0),
		_endRow(ms.nrow()),
		_threadCount(threadCount)
	{ }
	
	void InitializeInput(const class Model& model);
	
	DFTPredictionInput& Input() { return _dftInput; }
	
	~LMSPredicter();
	
	void Start();
	
	bool GetNextRow(RowData& data)
	{
		return _bufferedOutputLane.read(data);
	}
	
	void FinishRow(RowData& data)
	{
		_bufferedBufferLane.write(data);
	}
	
	std::mutex &IOMutex() { return _mutex; }
	
	void SetApplyBeam(bool applyBeam) { _applyBeam = applyBeam; }
	void SetStartRow(size_t startRow) { _startRow = startRow; }
	void SetEndRow(size_t endRow) { _endRow = endRow; }
	void SetUseModelColumn(bool useModelColumn) { _useModelColumn = useModelColumn; }
private:
	void ReadThreadFunc();
	void PredictThreadFunc();
	void clearBuffers();
	
	casacore::MeasurementSet &_ms;
	std::unique_ptr<LBeamEvaluator> _beamEvaluator;
	bool _applyBeam, _useModelColumn;
	
	DFTPredictionInput _dftInput;
	std::mutex _mutex;
	boost::barrier _barrier;
	boost::asio::io_service _ioService;
	
	const size_t _laneSize;
	ao::lane<RowData> _workLane, _outputLane, _availableBufferLane;
	lane_read_buffer<RowData> _bufferedOutputLane;
	lane_write_buffer<RowData> _bufferedBufferLane;
	
	std::unique_ptr<std::thread> _readThread;
	std::vector<std::thread> _workThreadGroup;
	std::unique_ptr<DFTPredictionAlgorithm> _predicter;
	std::vector<MC2x2*> _buffers;
	BandData _bandData;
	size_t _startRow, _endRow, _threadCount;
};

#endif
