// SPDX-License-Identifier: MIT
// NWKIT's external-process adapter. No IQ-TREE implementation is included here.
// A linked executable contains IQ-TREE and must be distributed under its GPL
// terms. NWKIT does not distribute that executable or the IQ-TREE library.

#include "tree/phylotree.h"
#include "model/modelfactory.h"
#include "model/modelsubst.h"
#include "utils/tools.h"
#include "utils/checkpoint.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef NWKIT_IQTREE_LIBRARY_SHA256
#error "Build this adapter with python -m nwkit.iqtree_library build"
#endif

namespace {
const char *prefix = "NWKIT_IQTREE_V1 ";

void leaf_names(Node *node, Node *parent, std::ostream &out) {
    if (node->isLeaf()) {
        if (node->name.empty() || node->name.find_first_of(" \t\r\n") != std::string::npos)
            throw std::runtime_error("Taxon names must be nonempty and whitespace-free");
        out << ' ' << node->name;
    }
    for (const auto edge : node->neighbors)
        if (edge->node != parent) leaf_names(edge->node, node, out);
}

void capabilities() {
    std::cout << "{\"protocol\":1,\"engine\":\"iqtree3\",\"library_version\":\""
              << iqtree_VERSION_MAJOR << '.' << iqtree_VERSION_MINOR << iqtree_VERSION_PATCH
              << "\",\"library_sha256\":\"" << NWKIT_IQTREE_LIBRARY_SHA256
              << "\",\"adapter_sha256\":\"" << NWKIT_IQTREE_ADAPTER_SHA256
              << "\",\"upstream_revision\":\"" << NWKIT_IQTREE_SOURCE_REVISION
              << "\"}" << std::endl;
}

void serve(PhyloTree &tree) {
    if (tree.rooted || !tree.getModel()->isReversible() || tree.getModel()->isMixture())
        throw std::runtime_error("A reversible, unrooted single-model tree is required");
    BranchVector branches;
    tree.getBranches(branches);
    std::cout << prefix << "READY " << branches.size() << std::endl;
    for (const auto edge : branches) {
        std::cout << prefix << "EDGE " << edge.first->findNeighbor(edge.second)->length;
        leaf_names(edge.first, edge.second, std::cout);
        std::cout << std::endl;
    }
    std::cout << prefix << "END" << std::endl;

    std::string line;
    while (std::getline(std::cin, line)) {
        if (line == "QUIT") return;
        std::istringstream request(line);
        std::string operation, extra;
        request >> operation;
        if (operation != "EVAL") throw std::runtime_error("Expected EVAL or QUIT");
        std::vector<double> lengths(branches.size());
        for (double &length : lengths)
            if (!(request >> length) || !std::isfinite(length) || length <= 0)
                throw std::runtime_error("One finite positive length per branch is required");
        if (request >> extra) throw std::runtime_error("Unexpected extra branch length");

        for (size_t i = 0; i < branches.size(); ++i) {
            const auto edge = branches[i];
            edge.first->findNeighbor(edge.second)->length = lengths[i];
            edge.second->findNeighbor(edge.first)->length = lengths[i];
        }
        // Invalidate branch-dependent partials using IQ-TREE's own interface.
        // The alignment, model and allocated likelihood buffers remain loaded.
        tree.clearAllPartialLH();
        tree.theta_computed = false;
        const double log_likelihood = tree.computeLikelihood();
        if (!std::isfinite(log_likelihood)) throw std::runtime_error("Nonfinite likelihood");
        std::vector<double> first(branches.size()), second(branches.size());
        for (size_t i = 0; i < branches.size(); ++i) {
            const auto edge = branches[i];
            tree.theta_computed = false;
            tree.computeLikelihoodDerv(
                static_cast<PhyloNeighbor *>(edge.first->findNeighbor(edge.second)),
                static_cast<PhyloNode *>(edge.first), &first[i], &second[i]);
            if (!std::isfinite(first[i]) || !std::isfinite(second[i]))
                throw std::runtime_error("Nonfinite branch derivative");
        }
        std::cout << prefix << "VALUE " << log_likelihood;
        for (double value : first) std::cout << ' ' << value;
        for (double value : second) std::cout << ' ' << value;
        std::cout << std::endl;
    }
}
}  // namespace

int main(int argc, char **argv) {
    if (argc == 2 && std::string(argv[1]) == "--capabilities") {
        capabilities();
        return 0;
    }
    try {
        auto &params = Params::getInstance();
        parseArg(argc, argv, params);
        // The portable SIMD kernel is available in both x86-64 and ARM64 builds.
        params.SSE = LK_SSE2;
        verbose_mode = VB_QUIET;
        if (!params.aln_file || !params.user_file || params.num_threads < 1)
            throw std::runtime_error("Alignment, fixed tree and positive thread count are required");
        init_random(params.ran_seed);
        Alignment alignment(params.aln_file, params.sequence_type, params.intype, params.model_name);
        Checkpoint checkpoint;
        PhyloTree tree(&alignment);
        tree.setParams(&params);
        tree.setCheckpoint(&checkpoint);
        bool rooted = false;
        tree.readTree(params.user_file, rooted);
        tree.setAlignment(&alignment);
        tree.setRootNode(params.root);
        std::unique_ptr<ModelsBlock> definitions(readModelsDefinition(params));
        tree.setModelFactory(new ModelFactory(params, params.model_name, &tree, definitions.get()));
        tree.setModel(tree.getModelFactory()->model);
        tree.setRate(tree.getModelFactory()->site_rate);
        tree.getModelFactory()->setCheckpoint(&checkpoint);
        tree.setLikelihoodKernel(params.SSE);
        tree.setNumThreads(params.num_threads);
        alignment.orderPatternByNumChars(PAT_VARIANT);
        tree.ensureNumberOfThreadsIsSet(nullptr);
        tree.initializeAllPartialLh();
        std::cout << std::setprecision(17);
        serve(tree);
    } catch (const std::exception &error) {
        std::cout << prefix << "ERROR " << error.what() << std::endl;
        return 1;
    } catch (const std::string &error) {
        std::cout << prefix << "ERROR " << error << std::endl;
        return 1;
    } catch (const char *error) {
        std::cout << prefix << "ERROR " << error << std::endl;
        return 1;
    }
    return 0;
}
