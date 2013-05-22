OBJECTST = General.o Layer.o Sequence.o Model.o NNt.o NN.o Ensemble.o Predict.o 
HEADERST = General.h Layer.h Model.h NNt.h NN.h Sequence.h Ensemble.h

CC	= g++
CFLAGS  = -O3

predict_seq_ss : $(OBJECTST) 
	$(CC)   $(CFLAGS) -o $@ $(OBJECTST) $(LLIBS)

Predict.o : Predict.cxx $(HEADERST)
	$(CC)  -c $(CFLAGS) -Wno-deprecated $*.cxx

Ensemble.o : Ensemble.cxx $(HEADERST)
	$(CC)   -c $(CFLAGS) -Wno-deprecated $*.cxx

General.o : General.cxx $(HEADERST)
	$(CC)   -c $(CFLAGS) -Wno-deprecated $*.cxx

Layer.o : Layer.cxx $(HEADERST)
	$(CC)   -c $(CFLAGS) -Wno-deprecated $*.cxx

Sequence.o : Sequence.cxx $(HEADERST)
	$(CC)   -c $(CFLAGS) -Wno-deprecated $*.cxx

Model.o : Model.cxx $(HEADERST)
	$(CC)   -c $(CFLAGS) -Wno-deprecated $*.cxx

NNt.o : NNt.cxx $(HEADERST)
	$(CC)   -c $(CFLAGS) -Wno-deprecated $*.cxx

NN.o : NN.cxx $(HEADERST)
	$(CC)   -c $(CFLAGS) -Wno-deprecated $*.cxx
	
clean:
	rm *.o predict 
