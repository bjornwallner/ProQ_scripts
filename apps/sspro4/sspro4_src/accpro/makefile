OBJECTST = General.o Layer.o Model.o NNt.o NN.o Ensemble.o PredictSeq.o 
HEADERST = General.h Layer.h Model.h NNt.h NN.h Sequence.h Ensemble.h

predict_seq_sa: $(OBJECTST) 
	g++ -O3 -static -Wno-deprecated -o $@ $(OBJECTST) $(LLIBS) 
PredictSeq.o: PredictSeq.cxx $(HEADERST)
	g++ -O3 -static -Wno-deprecated -c PredictSeq.cxx
Ensemble.o: Ensemble.cxx  $(HEADERST)
	g++ -O3 -static -Wno-deprecated -c Ensemble.cxx
General.o: General.cxx $(HEADERST) 
	g++  -O3 -static -Wno-deprecated -c General.cxx
Layer.o: Layer.cxx  $(HEADERST)
	g++  -O3 -static -Wno-deprecated -c Layer.cxx
NN.o: NN.cxx  $(HEADERST)
	g++  -O3 -static -Wno-deprecated -c NN.cxx
NNt.o: NNt.cxx  $(HEADERST)
	g++  -O3 -static -Wno-deprecated -c NNt.cxx
Model.o: Model.cxx  $(HEADERST)
	g++  -O3 -static -Wno-deprecated  -c Model.cxx

