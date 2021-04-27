#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
std::vector< std::string >  splitStrCpp(std::string x,std::string del) {
  
  std::vector< std::string > result;
  std::string s = x;
  std::size_t pos = 0;
  std::string token;
  while ((pos = s.find(del)) != std::string::npos) {
    token = s.substr(0, pos);
    result.push_back(token);
    s.erase(0, pos + del.length());
  }
  result.push_back(s);
  return result;
}





// [[Rcpp::export]]
bool innerConnection(std::string x,std::string y) {
  
  std::vector<std::string> x_split = splitStrCpp(x,"_");
  // Remove first and last clusters
  x_split.erase(x_split.begin());
  x_split.pop_back();
  
  std::vector<std::string> y_split = splitStrCpp(y,"_");
  
  // If cluster exists in the inner clusters, return true
  if ((std::find(x_split.begin(), x_split.end(),y_split[0])!=x_split.end()) |
      (std::find(x_split.begin(), x_split.end(),y_split[1])!=x_split.end())) {
    return true;
  } else  {
    return false;
  }
}


// [[Rcpp::export]]
bool cycle(std::string x,std::string y) {
  
  std::vector<std::string> x_split = splitStrCpp(x,"_");
  // Remove first and last clusters
  std::string x_first = x_split.front();
  std::string x_last = x_split.back();
  
  std::vector<std::string> y_split = splitStrCpp(y,"_");
  
  // Cycle is formed if the first or last element of y is the same 
  // as the first or last element of x
  if (((x_first == y_split[0]) & (x_last == y_split[1])) | 
      ((x_first == y_split[1]) & (x_last == y_split[0]))) {
    return true;
  } else  {
    return false;
  }
}

// [[Rcpp::export]]
bool connectable(std::string x,std::string y) {
  
  std::vector<std::string> x_split = splitStrCpp(x,"_");
  std::string x_first = x_split.front();
  std::string x_last = x_split.back();
  
  std::vector<std::string> y_split = splitStrCpp(y,"_");
  std::string y_first = y_split.front();
  std::string y_last = y_split.back();
  
  // x and y are connectable if either the first or last cluster of y 
  // are the same as the first or last cluster of x
  if ((x_first == y_first) | (x_first == y_last) | 
      (x_last == y_first) | (x_last == y_last)) {
    return true;
  } else  {
    return false;
  }
}

// [[Rcpp::export]]
bool replace(std::string& str, const std::string& from, const std::string& to) {
  size_t start_pos = str.find(from);
  if(start_pos == std::string::npos)
    return false;
  str.replace(start_pos, from.length(), to);
  return true;
}


// [[Rcpp::export]]
std::string connectGraphs(std::string x,std::string y) {
  
  std::vector<std::string> x_split = splitStrCpp(x,"_");
  std::string x_first = x_split.front();
  std::string x_last = x_split.back();
  
  std::vector<std::string> y_split = splitStrCpp(y,"_");
  std::string y_first = y_split.front();
  std::string y_last = y_split.back();
  
  std::string result;
  if (x_first == y_first) {
    std::vector<std::string> x_copy = x_split;
    std::vector<int> x_copy_int(x_copy.size());
    std::transform(x_copy.begin(), x_copy.end(), x_copy_int.begin(), [](const std::string& val)
    {
      return std::stoi(val);
    });
    std::reverse(x_copy_int.begin(), x_copy_int.end());
    std::vector<std::string> x_copy_string(x_copy.size());
    std::transform(x_copy_int.begin(), x_copy_int.end(), x_copy_string.begin(), [](const int& val)
    {
      return std::to_string(val);
    });
    std::string x_rev;
    for (std::vector<std::string>::const_iterator i = x_copy_string.begin(); i != x_copy_string.end(); ++i)
      x_rev += "_" + *i;
    x_rev.erase(0,1);
    
    result = x_rev + "_" + y;
    replace(result,"_"+x_first+"_","_");
  } else if (x_first == y_last) {
    result = y + "_" + x;
    replace(result,"_"+x_first+"_","_");
  } else if (x_last == y_first) {
    result = x + "_" + y;
    replace(result,"_"+x_last+"_","_");
  } else if (x_last == y_last) {
    std::vector<std::string> y_copy = y_split;
    std::vector<int> y_copy_int(y_copy.size());
    std::transform(y_copy.begin(), y_copy.end(), y_copy_int.begin(), [](const std::string& val)
    {
      return std::stoi(val);
    });
    std::reverse(y_copy_int.begin(), y_copy_int.end());
    std::vector<std::string> y_copy_string(y_copy.size());
    std::transform(y_copy_int.begin(), y_copy_int.end(), y_copy_string.begin(), [](const int& val)
    {
      return std::to_string(val);
    });
    std::string y_rev;
    for (std::vector<std::string>::const_iterator i = y_copy_string.begin(); i != y_copy_string.end(); ++i)
      y_rev += "_" + *i;
    y_rev.erase(0,1);
    
    
    result = x + "_" + y_rev;
    replace(result,"_"+x_last+"_","_");
  } else {
    result = "";
  }
  return result;
}





// [[Rcpp::export]]
std::string specialKruskal(NumericMatrix d) {
  
  std::vector<double> clusterDistances;
  std::vector<std::string> clusterPairs;
  
  int nr = d.nrow(), nc = d.ncol(), i = 0;
  
  if (nr == 2) {
    return "1_2";
  }
  
  
  for ( ; i < nr; i++) {
    int j = 0;
    for ( ; j < nc; j++) {
      if (j <= i) {
        continue;
      } else {
        std::string s = std::to_string(i+1) + "_" + std::to_string(j+1);
        clusterDistances.push_back(d(i,j));
        clusterPairs.push_back(s);
      }
    }
  }
  
  // Create pairs
  std::vector<std::pair<double, std::string>> d_pair_min;
  d_pair_min.reserve(clusterPairs.size());
  std::transform(clusterDistances.begin(), clusterDistances.end(), clusterPairs.begin(), std::back_inserter(d_pair_min),
                 [](double a, std::string b) { return std::make_pair(a, b); });
  
  // Sort into ascending order based on the distances
  sort(d_pair_min.begin(), d_pair_min.end()); 
  
  // Start the special case of Kruskal's algorithm
  std::vector<std::string> clusterGraph;
  for (size_t k = 0; k < d_pair_min.size() ; k++)
  {
    // std::cout << k << ' ' << std::endl;
    std::pair<int, std::string> min_pair_k = d_pair_min[k];
    std::string pair_min = min_pair_k.second;
    
    // std::cout << pair_min << ' ' << std::endl;
    
    
    // If the graph is empty, initialize it with the pair 
    // having the shortest distance
    if (clusterGraph.empty()) {
      clusterGraph.push_back(pair_min);
      continue;
    }

    // Check if the graph is ready: all clusters are connected together
    
    if (clusterGraph.size() == 1) {
      std::vector < std::string > clusterGraph_1_split = splitStrCpp(clusterGraph[0],"_");
      size_t clusters = d.ncol();
      if (clusterGraph_1_split.size() == clusters)
      {
        return clusterGraph[0];
      }
    }
    
    // Check if pair_min can be linked to any of
    // the existing disconnected graph(s)
    // 1. Inner connections are forbidden.
    // 2. Connections that make a disconnected subgraph cycle are forbidden.
    // 3. Else, the new link will be added to an existing subgraph.
    // 4. If the connection can not be linked to any existing disconnected
    //    subgraph, then the new connection will form a new disconnected subgraph.
    
    bool illegal = false;
    for (size_t i = 0; i < clusterGraph.size() ; i++) {
      // Test conditions 1 and 2
      if (innerConnection(clusterGraph[i],pair_min) |
          cycle(clusterGraph[i],pair_min)) {
        illegal = true;
        break;
      }
    }
    
    
    if (illegal)
    {
      continue;
    }
    std::string connectedSubGraph;
    // Test condition 3 and link to existing subgraph
    bool notConnectable = true;
    for (size_t i = 0; i < clusterGraph.size() ; i++) {
      if (connectable(clusterGraph[i],pair_min))
      {
        connectedSubGraph = connectGraphs(clusterGraph[i],pair_min);
        clusterGraph[i] = connectedSubGraph;
        notConnectable = false;
        break;
      }
    }
    

    // If a new connection was made, then
    // check if after adding the new connection some of the subgraphs
    // become connectable. If yes, connect them.
    
    if (!notConnectable)
    {

      for (size_t i = 0; i < clusterGraph.size() ; i++) {
        std::string clusterGraph_i = clusterGraph[i];
        if (connectedSubGraph != clusterGraph_i)
        {
          if (connectable(clusterGraph_i,connectedSubGraph)) {

            std::string newConnectedGraph = connectGraphs(clusterGraph_i,connectedSubGraph);
            clusterGraph.push_back(newConnectedGraph);
            clusterGraph.erase(std::remove(clusterGraph.begin(), clusterGraph.end(), clusterGraph_i), clusterGraph.end());
            clusterGraph.erase(std::remove(clusterGraph.begin(), clusterGraph.end(), connectedSubGraph), clusterGraph.end());
            break;
          }
        }
      }
    } else {
      // 4: Add as a new disconnected subgraph
      clusterGraph.push_back(pair_min);
    }
    
  }
  return "";
}
