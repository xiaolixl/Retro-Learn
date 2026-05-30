from pathlib import Path
from typing import Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from agent_runtime import RetrosynthesisAgent

PROJECT_ROOT = Path(__file__).resolve().parent
STATIC_DIR = PROJECT_ROOT / "static"
USER_OUTPUT_DIR = PROJECT_ROOT / "user_output"

app = FastAPI(
    title="SimpRetro LLM Agent API",
    description="Natural-language retrosynthesis agent built on top of the SimpRetro backend.",
    version="0.2.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class AgentQueryRequest(BaseModel):
    message: str = Field(..., description="Natural-language user request")
    model: Optional[str] = Field(default=None, description="Optional OpenAI model override")


@app.get("/health")
def healthcheck() -> dict:
    return {"status": "ok"}


@app.post("/agent/query")
def query_agent(request: AgentQueryRequest) -> dict:
    try:
        agent = RetrosynthesisAgent(model=request.model)
        result = agent.run(request.message)

        if result.get("status") == "success" and "output_dir" in result:
            output_dir = Path(result["output_dir"])
            relative = output_dir.relative_to(USER_OUTPUT_DIR).as_posix()
            result["image_base_url"] = f"/user_output/{relative}"
            if result.get("viz_html"):
                viz_path = Path(result["viz_html"])
                if viz_path.exists():
                    viz_rel = viz_path.relative_to(USER_OUTPUT_DIR).as_posix()
                    result["viz_html_url"] = f"/user_output/{viz_rel}"

        return result
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc


# Serve the SPA frontend at /, and user output images at /user_output.
if USER_OUTPUT_DIR.exists():
    app.mount("/user_output", StaticFiles(directory=str(USER_OUTPUT_DIR)), name="user_output")

if STATIC_DIR.exists():
    from fastapi.responses import FileResponse

    @app.get("/")
    def serve_frontend():
        return FileResponse(str(STATIC_DIR / "index.html"))

    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")
