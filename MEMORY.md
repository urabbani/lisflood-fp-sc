# MEMORY.md - Long-term Memory for SimOne

Last updated: 2026-03-28 (19:40)

## User Profile

**Name:** Dr. Umair Rabbani  
**Contact:** +923354110084  
**Timezone:** Asia/Karachi (GMT+5)  
**Workspace:** /home/umair/.openclaw/workspace (WSL2 on Windows)

## Infrastructure

**Machine:** Omen (192.168.4.165)  
**OS:** Linux 6.6.87.2-microsoft-standard-WSL2 (x64)  
**Node:** v23.11.1

### Ollama Server
- **Host:** 10.0.0.205:11434
- **GPU:** RTX 5090 with 24GB VRAM
- **Available models:**
  - qwen3.5:35b-a3b-q8_0 - Main reasoning model (262k context, vision-enabled)
  - qwen3-embedding:latest - 4096-dim embeddings for semantic memory

### OpenClaw Setup
- **Version:** 2026.3.24 (installed version)
- **Gateway:** Running on ws://127.0.0.1:18789
- **WhatsApp Number:** +923174160617 (SimOne's number)
- **Connected to:** +923354110084 (Dr. Umair Rabbani)
- **Systemd:** Enabled and active

## Memory System

**Plugin:** openclaw-memory-ollama (community plugin)  
**Database:** LanceDB at ~/.openclaw/memory/lancedb  
**Embedding:** qwen3-embedding:latest (4096 dimensions)  
**Features:**
- Semantic search (vector similarity)
- Auto-capture enabled (scans for facts, preferences, decisions)
- Auto-recall enabled (injects relevant memories into context)
- 100% local, zero cloud dependencies

**Memory files:** 5 daily files (2026-03-24 to 2026-03-28)

## Projects & Work

### MolmoWeb (2026-03-26)
- Downloaded locally in workspace/molmoweb/
- Download paused (~4GB of 8GB)
- **Decision:** Will build Ollama-based web agent instead
- **Reasoning:** Ollama has RTX 5090, saves 4GB, faster inference
- qwen3.5:35b-a3b-q8_0 confirmed to support image processing

### BMS Battery System Analysis (2026-03-27)
**Battery Pack:** 15S LiFePO4 (15 cells in series)
- Full capacity: 280.0Ah
- Nominal voltage range: 3.2-3.65V per cell
- Pack voltage: 55.5-57.0V
- Temperature range: 25-33°C
- BMS: JK Tech app

**Critical Issues Identified:**

1. **Cell 16 Imbalance (SEVERE)**
   - Voltage gap: Cell 16 ~3.70V vs lowest cell ~3.32V (0.38V difference)
   - Overcharge protection events: 2,000+ occurrences in 16 days
   - Pattern: Cell 16 overcharges → Protection ON → Charging OFF → Voltage drops → Protection OFF → Cycle repeats every 20-40 seconds
   - Impact: Prevents full pack charging, reduces available capacity

2. **Coprocessor Communication Errors (RECURRING)**
   - Occurrences: 4 events on 2026-03-14, 03-16, 03-17, 03-23
   - Symptoms:
     * All temperature sensors jump to 56°C simultaneously
     * Battery voltage reading drops
     * Charging turns OFF automatically
     * Error: "Abnormal coprocessor communication"
   - Likely causes: BMS master module issue, balancing board malfunction, wiring/connector problem

**Immediate Actions Recommended:**
1. Check Cell 16 physical connection (loose terminal, poor contact)
2. Measure individual cell voltage with multimeter
3. Inspect for cell swelling or damage
4. Check BMS board connectors and wiring
5. Consider external balancer for severe imbalance
6. BMS firmware update if available

### Annotations Feature Planning (2026-03-27)
**Status:** Planning phase - detailed implementation documents created
**Goal:** Collaborative drawing & annotation system (Google Earth-style)
**Key Features:**
- Draw points, lines, polygons on map
- Add title, description, priority, status, tags
- Save to PostgreSQL/PostGIS with full attribution (created by, timestamp)
- Export to GeoJSON/KML
- Collaborative: All users see each other's annotations

**Files Created:**
- `user-annotations-feature-plan.md` - Full implementation plan (42KB)
- `github-issue-annotations-feature.md` - Detailed GitHub issue for Claude Code (19.5KB)
- `create_annotations_table.sql` - Database migration script
- `annotations.mjs` - Complete Express API server (port 3002)
- `floodrisk-annotations-api.service` - Systemd service template
- `annotations-quick-start.md` - Step-by-step setup guide

**Architecture:**
- Database: `public.user_annotations` table in PostgreSQL/PostGIS
- Backend: Express API on port 3002 (separate from Impact API on 3001)
- Frontend: OpenLayers vector layer, React components (Toolbar, Editor, List)
- Export: GeoJSON and KML formats
- Styling: Priority-based color coding (green/yellow/orange/red)

**Timeline:** 17-23 hours (2-3 days) for full implementation

### WhatsApp Integration (2026-03-24/27)
- Gateway initially had cycling issues (status 499) - 2026-03-25
- Resolved with restart - 2026-03-25
- **Major Crisis (2026-03-27):** Gateway cycling every 60 seconds - RECURRING ISSUE
  - Root cause: Corrupted credentials file (creds.json)
  - Pattern: Cycles for hours → stabilizes 20-50 min → cycles again
  - Timeline:
    * Cycling 22:06-05:02 (7+ hours) → stopped at 08:13
    * Cycling resumed 09:03-09:13 (10 min)
    * Reconnected at 09:48:42 - currently stable (~25 min)
  - Temporary fixes: Gateway self-recovers but issue returns
  - **Permanent fix needed:** Delete corrupted credentials and re-link
  - Commands: `rm ~/.openclaw/credentials/whatsapp/default/creds.json*` then `openclaw channels login --channel whatsapp --account default`
  - Credentials path: ~/.openclaw/credentials/whatsapp/default/creds.json

- **Chronic Pattern (2026-03-28):** Credentials corruption occurring predictably every 35 minutes
  - Pattern consistent throughout day: 06:44, 07:19, 07:54, 08:29, 09:04, 09:39 (35 min intervals)
  - Continued: 10:15, 10:50, 11:25, 12:00, 12:35 (35 min intervals)
  - Brief anomaly at 12:48 (13 min) then returned to 35 min pattern
  - Auto-restore from backup working flawlessly - gateway stays connected
  - No message delivery issues
  - Root cause: Underlying WhatsApp Web session instability
  - Status: MONITORING - Auto-restore reliable, predictable pattern managed
  - Recommended: Reinstall WhatsApp channel when convenient for permanent fix
- All messages delivering correctly
- WhatsApp linked: +923174160617

## Preferences

- **Local-first:** Prefers local models over cloud APIs
- **Privacy-focused:** Wants 100% local memory (no cloud embeddings)
- **GPU-powered:** Utilizing RTX 5090 on Ollama server
- **WSL2:** Running on Windows Subsystem for Linux

## Technical Notes

### OpenClaw Plugins
- ✅ lossless-claw (LCM context management)
- ✅ memory-ollama (Ollama-based semantic memory)
- ✅ whatsapp (messaging channel)
- ✅ memory-core (disabled, replaced by memory-ollama)
- ✅ memory-lancedb (disabled, replaced by memory-ollama)

### Issues Resolved
1. WhatsApp reply routing (2026-03-24)
2. Gateway cycling (2026-03-25)
3. Memory indexing - switched to Ollama embeddings (2026-03-26)
4. FTS5 unavailable in Node runtime (using LIKE fallback)
5. **WhatsApp credential corruption crisis (2026-03-27-28):**
   - Gateway cycling every 60 seconds for 7+ hours (2026-03-27)
   - Root cause: Corrupted creds.json file
   - Solution: Gateway restart cleared corruption
   - Troubleshooting key: Check logs for "restored corrupted WhatsApp creds.json from backup"
   - Evolved to chronic pattern (2026-03-28): Credentials corrupt predictably every ~35 minutes
   - Auto-restore from backup working reliably
   - Status: MONITORING - Pattern predictable, managed, no message delivery issues
   - Recommended: Reinstall WhatsApp channel when convenient for permanent fix

### Known Issues
- FTS5 module unavailable in current Node runtime (full-text search falls back to LIKE)
- memory-core required cloud embeddings (solved with memory-ollama plugin)
- **WhatsApp credentials corruption (2026-03-28):** Chronic predictable pattern every ~35 minutes; auto-restore working reliably; recommended to reinstall WhatsApp channel for permanent fix

## Important Decisions

1. **Ollama over Cloud:** Using local models for both reasoning and embeddings
2. **Memory Privacy:** Chose memory-ollama to keep all data local
3. **Web Browsing:** Will build custom agent using qwen3.5 + Playwright instead of MolmoWeb
4. **GPU Utilization:** Prioritizing RTX 5090 for heavy workloads

---

*This file is maintained by SimOne and updated as new information is learned.*
