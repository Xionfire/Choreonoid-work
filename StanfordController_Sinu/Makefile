CXXFLAGS += -fPIC -std=c++14 `pkg-config --cflags choreonoid`
CONTROLLER = StanfordController_sinu.so

StanfordController_sinu.so: StanfordController_sinu.o
	g++ -shared -o StanfordController_sinu.so StanfordController_sinu.o `pkg-config --libs choreonoid` -lCnoidBody
