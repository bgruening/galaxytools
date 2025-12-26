from cwl_utils.parser import cwl_v1_0, cwl_v1_1, cwl_v1_2, load_document_by_uri

from .util import get_path


def test_cuda_requirement_v1_0() -> None:
    """Test that CUDARequirement objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/cuda-requirement_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_0.CUDARequirement)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:CUDARequirement"


def test_cuda_requirement_v1_1() -> None:
    """Test that CUDARequirement objects are correctly loaded for CWL v1.1."""
    uri = get_path("testdata/extensions/cuda-requirement_v1_1.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_1.CUDARequirement)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:CUDARequirement"


def test_cuda_requirement_v1_2() -> None:
    """Test that CUDARequirement objects are correctly loaded for CWL v1.2."""
    uri = get_path("testdata/extensions/cuda-requirement_v1_2.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_2.CUDARequirement)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:CUDARequirement"


def test_load_listing_requirement_v1_0() -> None:
    """Test that LoadListingRequirement objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/load-listing-requirement_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_0.LoadListingRequirement)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:LoadListingRequirement"


def test_loop_v1_2() -> None:
    """Test that Loop and LoopInput objects are correctly loaded for CWL v1.2."""
    uri = get_path("testdata/extensions/single-var-loop_v1_2.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    cwl_step = next(iter(cwl_obj.steps))
    loop_req = next(iter(cwl_step.requirements))
    assert isinstance(loop_req, cwl_v1_2.Loop)
    assert isinstance(next(iter(loop_req.loop)), cwl_v1_2.LoopInput)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["steps"][0]["requirements"][0]["class"] == "cwltool:Loop"


def test_inplace_update_requirement_v1_0() -> None:
    """Test that InplaceUpdateRequirement objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/inplace-update-requirement_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(
        next(iter(cwl_obj.requirements)), cwl_v1_0.InplaceUpdateRequirement
    )
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:InplaceUpdateRequirement"


def test_mpi_requirement_v1_0() -> None:
    """Test that MPIRequirement objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/mpi-requirement_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_0.MPIRequirement)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:MPIRequirement"


def test_mpi_requirement_v1_1() -> None:
    """Test that MPIRequirement objects are correctly loaded for CWL v1.1."""
    uri = get_path("testdata/extensions/mpi-requirement_v1_1.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_1.MPIRequirement)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:MPIRequirement"


def test_mpi_requirement_v1_2() -> None:
    """Test that MPIRequirement objects are correctly loaded for CWL v1.2."""
    uri = get_path("testdata/extensions/mpi-requirement_v1_2.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_2.MPIRequirement)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:MPIRequirement"


def test_network_access_v1_0() -> None:
    """Test that NetworkAccess objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/network-access_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_0.NetworkAccess)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:NetworkAccess"


def test_process_generator_v1_0() -> None:
    """Test that ProcessGenerator objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/process-generator_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(cwl_obj, cwl_v1_0.ProcessGenerator)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["class"] == "cwltool:ProcessGenerator"


def test_process_generator_v1_1() -> None:
    """Test that ProcessGenerator objects are correctly loaded for CWL v1.1."""
    uri = get_path("testdata/extensions/process-generator_v1_1.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(cwl_obj, cwl_v1_1.ProcessGenerator)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["class"] == "cwltool:ProcessGenerator"


def test_process_generator_v1_2() -> None:
    """Test that ProcessGenerator objects are correctly loaded for CWL v1.2."""
    uri = get_path("testdata/extensions/process-generator_v1_2.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(cwl_obj, cwl_v1_2.ProcessGenerator)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["class"] == "cwltool:ProcessGenerator"


def test_secrets_v1_0() -> None:
    """Test that Secrets objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/secrets_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_0.Secrets)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:Secrets"


def test_secrets_v1_1() -> None:
    """Test that Secrets objects are correctly loaded for CWL v1.1."""
    uri = get_path("testdata/extensions/secrets_v1_1.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_1.Secrets)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:Secrets"


def test_secrets_v1_2() -> None:
    """Test that Secrets objects are correctly loaded for CWL v1.2."""
    uri = get_path("testdata/extensions/secrets_v1_2.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_2.Secrets)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:Secrets"


def test_shm_size_v1_0() -> None:
    """Test that ShmSize objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/shm-size_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_0.ShmSize)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:ShmSize"


def test_shm_size_v1_1() -> None:
    """Test that ShmSize objects are correctly loaded for CWL v1.1."""
    uri = get_path("testdata/extensions/shm-size_v1_1.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_1.ShmSize)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:ShmSize"


def test_shm_size_v1_2() -> None:
    """Test that ShmSize objects are correctly loaded for CWL v1.2."""
    uri = get_path("testdata/extensions/shm-size_v1_2.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_2.ShmSize)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:ShmSize"


def test_time_limit_v1_0() -> None:
    """Test that TimeLimit objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/time-limit_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_0.TimeLimit)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:TimeLimit"


def test_work_reuse_v1_0() -> None:
    """Test that WorkReuse objects are correctly loaded for CWL v1.0."""
    uri = get_path("testdata/extensions/work-reuse_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert isinstance(next(iter(cwl_obj.requirements)), cwl_v1_0.WorkReuse)
    cwl_dict = cwl_obj.save(top=True)
    assert cwl_dict["requirements"][0]["class"] == "cwltool:WorkReuse"
