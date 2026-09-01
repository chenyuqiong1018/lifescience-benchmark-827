from pathlib import Path
import json
import shutil
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm


def main(input_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    raw = pd.read_excel(input_dir / "data.xlsx", sheet_name="Sheet1")
    cols = ["Efficacy", "Age", "Gender", "BMI"]
    data = raw[cols].copy()
    data["Age"] = pd.to_numeric(data["Age"], errors="coerce")
    data["BMI"] = pd.to_numeric(data["BMI"], errors="coerce")
    data = data.dropna(subset=cols)
    data = data[data["Efficacy"].isin(["PR", "SD", "PD"]) & data["Gender"].isin(["Female", "Male"])].copy()
    data["response"] = data["Efficacy"].eq("PR").astype(int)
    data["Gender_Male"] = data["Gender"].eq("Male").astype(int)
    X = sm.add_constant(data[["BMI", "Age", "Gender_Male"]].astype(float), has_constant="add")
    model = sm.Logit(data["response"], X).fit(disp=False)
    coef = pd.DataFrame({
        "term": ["Intercept", "BMI", "Age", "Gender[Male]"],
        "estimate": model.params.to_numpy(),
        "std_error": model.bse.to_numpy(),
        "z": model.tvalues.to_numpy(),
        "p_value": model.pvalues.to_numpy(),
        "odds_ratio": np.exp(model.params.to_numpy()),
    })
    coef.to_csv(output_dir / "model_coefficients.csv", index=False)
    metadata = {
        "outcome": "Efficacy",
        "outcome_coding": {"PR": 1, "SD": 0, "PD": 0},
        "gender_reference_level": "Female",
        "gender_contrast": "Male vs Female",
        "predictors": ["BMI", "Age", "Gender"],
        "complete_case_variables": cols,
        "n_input_rows": int(len(raw)),
        "n_complete_cases": int(len(data)),
        "n_response_1": int(data.response.sum()),
        "n_response_0": int((1 - data.response).sum()),
        "model": "binary logistic regression (maximum likelihood)",
    }
    (output_dir / "model_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    shutil.copy2(Path(__file__), output_dir / "analysis.py")
    age = coef.loc[coef.term.eq("Age")].iloc[0]
    report = f"""# Treatment-response logistic regression

Efficacy was coded as response (`PR=1`) versus no response (`SD/PD=0`). Complete cases were selected only on Efficacy, Age, Gender, and BMI. Female is the gender reference level, so `Gender[Male]` compares Male with Female. The fitted maximum-likelihood logistic regression used BMI, age, and gender.

There were {len(data)} complete cases ({int(data.response.sum())} responses and {int((1-data.response).sum())} non-responses). The age log-odds coefficient was {age.estimate:.6f} (SE {age.std_error:.6f}, z {age.z:.6f}, two-sided p={age.p_value:.6g}; OR={age.odds_ratio:.6f} per year).
"""
    (output_dir / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
