[x , fs ] = jukebox;

%%
pwelch (x , [] , [] , [] , fs );
%%
sound (x , fs);
%%
freq = 1000*[.1 1.1 2.45];
freq = [-fliplr(freq) freq];
r = .1;
z = exp ( -1j*2 * pi * freq / fs);
p = exp ( -1j*2 * pi * freq / fs ) * r;
b = poly ( z );
a = poly ( p );
%%
y = filter (b , a , x);
pwelch (y , [] , [] , [] , fs);
%%
sound (y , fs);