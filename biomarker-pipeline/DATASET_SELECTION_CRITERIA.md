# Dataset Selection Criteria for CAMK2D Meta-Analysis

This document defines the systematic criteria for selecting and validating Gene Expression Omnibus (GEO) datasets for CAMK2D biomarker research and meta-analysis.

## Overview

The CAMK2D biomarker pipeline implements a rigorous, multi-tier validation system to ensure only high-quality, relevant datasets are included in meta-analysis. This approach maximizes statistical power while maintaining scientific rigor.

## Primary Selection Criteria

### 1. **CAMK Gene Coverage** (Weight: 40%)
- **Critical Requirement**: Dataset must contain measurable expression data for CAMK genes
- **Primary Target**: CAMK2D must be detectable (10 points)
- **Secondary Targets**: Other CAMK family genes (CAMK2A, CAMK2B, CAMK2G, CAMK1, CAMK4)
- **Minimum Threshold**: At least 1 CAMK gene detectable
- **Scoring**: 
  - CAMK2D detected: 10/10 points
  - Other CAMK genes only: 5-8/10 points
  - No CAMK genes: 0/10 points (exclusion)

### 2. **Dataset Quality** (Weight: 30%)
- **Sample Size**: Minimum 10 samples, optimal ≥50 samples
- **Gene Coverage**: Minimum 5,000 genes, optimal ≥20,000 genes
- **Missing Data**: Maximum 30% missing values
- **Platform Quality**: RNA-seq preferred over microarray
- **Batch Effects**: Assessment and correction capability
- **Scoring Formula**: (Sample_Score + Coverage_Score + Missing_Penalty) / 3

### 3. **Research Relevance** (Weight: 30%)
- **Cardiac Relevance** (30% of relevance score)
  - Heart failure, atrial fibrillation, cardiac hypertrophy: 10/10 points
  - Other cardiac conditions: 7-9/10 points
  - Non-cardiac but cardiovascular: 5-6/10 points
  - Non-cardiovascular: 0-4/10 points

- **CAMK2D Focus** (30% of relevance score)
  - Direct CAMK2D manipulation/study: 10/10 points
  - Calcium signaling focus: 7-9/10 points
  - General cardiac study: 5-6/10 points
  - Unrelated research: 0-4/10 points

- **Species Relevance** (20% of relevance score)
  - Human or Mouse: 10/10 points
  - Other mammalian models: 6-8/10 points
  - Non-mammalian: 3-5/10 points

- **Publication Quality** (20% of relevance score)
  - Based on journal impact, study design, sample size, methodology

## Dataset Categories and Priorities

### Tier 1: CAMK2D-Specific Studies (Highest Priority)
- **Examples**: GSE225336 (UK Biobank fibrosis), GSE266652 (gene editing), GSE201018 (splicing)
- **Characteristics**: Direct CAMK2D focus, cardiac tissue, high sample quality
- **Requirements**: Overall score ≥8.0, CAMK2D detection required
- **Usage**: Core datasets for primary analysis

### Tier 2: Cardiac CAMK Studies (High Priority) 
- **Examples**: Heart failure, atrial fibrillation with CAMK expression
- **Characteristics**: Cardiac focus with CAMK relevance
- **Requirements**: Overall score ≥7.0, cardiac tissue, CAMK detection
- **Usage**: Supporting datasets for meta-analysis

### Tier 3: Cardiac General Studies (Moderate Priority)
- **Examples**: General cardiac disease studies with CAMK expression
- **Characteristics**: Cardiac tissue but not CAMK-focused
- **Requirements**: Overall score ≥6.0, cardiac tissue
- **Usage**: Validation and context datasets

### Tier 4: Cross-Species Validation (Species-Specific Priority)
- **Examples**: Mouse cardiac models with human ortholog mapping
- **Characteristics**: Model organisms with translational relevance
- **Requirements**: Overall score ≥6.0, ortholog mapping available
- **Usage**: Cross-species validation studies

## Validation Workflow

### Step 1: Initial Screening
```r
# Check basic eligibility
- GEO dataset accessible
- Species: Human or Mouse preferred
- Tissue: Cardiac preferred
- Sample size: ≥10 samples
```

### Step 2: CAMK Detection Validation
```r
# Verify CAMK gene presence
validate_camk_detection(gse_id, expression_data)
- CAMK2D presence: Critical
- Other CAMK genes: Beneficial
- Gene symbol mapping: Verified
```

### Step 3: Quality Assessment
```r
# Evaluate dataset quality
assess_dataset_quality(gse_id, expression_data, metadata)
- Sample size adequacy
- Gene coverage sufficiency  
- Missing data assessment
- Platform compatibility
```

### Step 4: Relevance Scoring
```r
# Score research relevance
score_research_relevance(gse_id)
- Cardiac disease relevance
- CAMK2D research focus
- Species appropriateness
- Publication quality
```

### Step 5: Final Validation
```r
# Overall validation decision
calculate_validation_score(camk_results, quality_results, relevance_results)
- Weighted scoring algorithm
- Threshold-based inclusion
- Recommendation generation
```

## Inclusion/Exclusion Decision Tree

```
Dataset → CAMK Gene Detection?
    ├─ No → EXCLUDE (Critical requirement)
    └─ Yes → Quality Score ≥6.0?
        ├─ No → EXCLUDE (Quality threshold)
        └─ Yes → Overall Score ≥6.0?
            ├─ No → EXCLUDE (Combined threshold)
            └─ Yes → Relevance Assessment
                ├─ Score ≥8.0 → TIER 1 (Primary analysis)
                ├─ Score 7.0-7.9 → TIER 2 (Meta-analysis)
                ├─ Score 6.0-6.9 → TIER 3 (Validation)
                └─ Score <6.0 → EXCLUDE
```

## Special Considerations

### Recent High-Value Studies
- **GSE225336**: UK Biobank (41,505 participants) - Population genetics
- **GSE266652**: CAMK2D base editing - Therapeutic relevance
- **GSE201018**: CAMK2D splicing variants - Mechanistic insights
- **GSE228762**: CAMK2D CRISPR knockout - Functional validation

### Platform Considerations
- **RNA-seq**: Preferred for comprehensive expression profiling
- **Microarray**: Acceptable for focused gene analysis
- **Amplicon-seq**: Specialized applications (gene editing validation)
- **Single-cell**: Future consideration for cellular heterogeneity

### Cross-Species Mapping
- Human-Mouse ortholog validation required
- Gene symbol harmonization critical
- Platform-specific probe mapping
- Evolutionary conservation assessment

## Quality Control Checkpoints

### Pre-Analysis Validation
1. **Dataset Accessibility**: GEO download successful
2. **Metadata Completeness**: Sample annotations available
3. **Expression Matrix Integrity**: No corruption, proper format
4. **Gene Symbol Mapping**: Successful symbol conversion

### Post-Loading Validation
1. **CAMK Gene Expression**: Detectable expression levels
2. **Sample Quality**: Outlier detection and removal
3. **Batch Effect Assessment**: Platform and study effects
4. **Missing Data Handling**: Imputation strategies

### Meta-Analysis Validation
1. **Cross-Study Harmonization**: Batch correction success
2. **Effect Size Homogeneity**: Heterogeneity assessment (I²)
3. **Publication Bias Testing**: Funnel plots, Egger's test
4. **Sensitivity Analysis**: Robustness to dataset exclusion

## Implementation Tools

### Validation Functions
- `validate_dataset()`: Comprehensive dataset validation
- `batch_validate_datasets()`: Multiple dataset validation
- `search_datasets()`: Criteria-based dataset search
- `get_dataset_stats()`: Configuration statistics

### Configuration Management
- `add_dataset()`: Add new datasets to configuration
- `get_priority_datasets()`: Retrieve high-priority datasets
- `get_all_dataset_configs()`: Access all configurations
- `get_analysis_config()`: Analysis parameters

## Continuous Improvement

### Regular Review Process
- **Monthly**: New dataset screening and addition
- **Quarterly**: Validation criteria refinement
- **Annually**: Comprehensive criteria review and update

### Feedback Integration
- User feedback incorporation
- Literature review updates
- Technology advancement adaptation
- Statistical method improvements

## References and Standards

### Methodological Standards
- PRISMA guidelines for systematic reviews
- MIAME standards for microarray experiments  
- MINSEQE standards for RNA-seq experiments
- FAIR data principles compliance

### Statistical Frameworks
- Random effects meta-analysis models
- Heterogeneity assessment (I², Q-test)
- Publication bias detection methods
- Cross-species validation approaches

---

*This document is maintained as part of the CAMK2D biomarker pipeline and updated regularly to reflect best practices and new developments in the field.*

**Version**: 1.0  
**Last Updated**: `r Sys.Date()`  
**Next Review**: Quarterly