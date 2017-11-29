APPS = adjmat2wpairs cluster_communities getting_communities compute_tanimoto
all: $(APPS)

adjmat2wpairs: adjmat2wpairs.cpp
	g++ -o $@ $<

cluster_communities: cluster_communities.cpp
	g++ -o $@ $<

getting_communities: getting_communities.cpp
	g++ -o $@ $<

compute_tanimoto: compute_tanimoto.cpp
	g++ -o $@ $<

clean:
	rm -rf $(APPS)

