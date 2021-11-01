CXXFLAGS += -fPIC -std=c++14 `pkg-config --cflags choreonoid`
CONTROLLER = StanfordController_3D.so

StanfordController_3D.so: StanfordController_3D.o
	g++ -shared -o StanfordController_3D.so StanfordController_3D.o `pkg-config --libs choreonoid` -lCnoidBody
