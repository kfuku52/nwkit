import builtins
import io
import os
import sqlite3
import stat
import struct
import tarfile
import zlib

import pytest
import requests

from nwkit import image as image_module
from nwkit.image import (
    MediaDownloadError,
    download_media,
    image_main,
    postprocess_media_file,
)
from tests.image_test_support import (
    DummyProvider,
    DummySession,
    make_image_args,
    read_tsv,
    write_valid_test_media,
)


class TestDownloadMedia:
    def test_download_media_uses_media_accept_header(self, tmp_path):
        destination = tmp_path / 'out' / 'image.bin'
        request_kwargs = {}

        class JPEGResponse:
            headers = {'Content-Type': 'image/jpeg'}
            url = 'https://example.org/download'

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'\xff\xd8\xff\xe0fakejpeg'

            def close(self):
                return None

        class JPEGSession:
            def get(self, *args, **kwargs):
                request_kwargs.update(kwargs)
                return JPEGResponse()

        download_media(
            JPEGSession(),
            'https://example.org/download',
            str(destination),
        )

        assert request_kwargs['headers']['Accept'].startswith('image/*')
        assert request_kwargs['headers']['Accept'] != 'application/json'
        assert image_module.HTTP_HEADERS is image_module.API_HTTP_HEADERS
        assert image_module.API_HTTP_HEADERS['Accept'] == 'application/json'

    def test_download_media_prefers_shared_cache_when_present(self, tmp_path):
        destination = tmp_path / 'out' / 'image.svg'
        cache_path = tmp_path / 'cache' / 'image.svg'
        cache_path.parent.mkdir(parents=True)
        cache_path.write_bytes(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')

        class FailingSession:
            def get(self, *args, **kwargs):
                raise AssertionError('network should not be used when cache exists')

        result = download_media(FailingSession(), 'https://images.example.org/item.svg', str(destination), cache_path=str(cache_path))

        assert result['status'] == 'cached'
        assert destination.read_bytes().startswith(b'<svg')

    def test_cached_destination_uses_umask_and_preserves_existing_mode(
        self,
        tmp_path,
    ):
        destination = tmp_path / 'out' / 'image.svg'
        cache_path = tmp_path / 'cache' / 'image.svg'
        cache_path.parent.mkdir(parents=True)
        cache_path.write_bytes(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')

        class FailingSession:
            def get(self, *args, **kwargs):
                raise AssertionError('network should not be used when cache exists')

        previous_umask = os.umask(0o027)
        try:
            download_media(
                FailingSession(),
                'https://images.example.org/item.svg',
                str(destination),
                cache_path=str(cache_path),
            )
        finally:
            os.umask(previous_umask)
        if os.name != 'nt':
            assert stat.S_IMODE(destination.stat().st_mode) == 0o640

        destination.chmod(0o664)
        download_media(
            FailingSession(),
            'https://images.example.org/item.svg',
            str(destination),
            cache_path=str(cache_path),
            reuse_destination=False,
        )
        if os.name != 'nt':
            assert stat.S_IMODE(destination.stat().st_mode) == 0o664

    def test_download_media_removes_partial_temp_file_on_error(self, tmp_path):
        destination = tmp_path / 'out' / 'image.svg'
        destination.parent.mkdir(parents=True)

        class FailingSession:
            def get(self, *args, **kwargs):
                raise requests.ConnectionError('connection dropped')

        with pytest.raises(requests.ConnectionError, match='connection dropped'):
            download_media(FailingSession(), 'https://images.example.org/item.svg', str(destination))

        assert not destination.exists()
        assert not (tmp_path / 'out' / 'image.svg.tmp').exists()

    def test_download_media_uses_content_type_to_resolve_extension(self, tmp_path):
        destination = tmp_path / 'out' / 'image.bin'

        class JPEGResponse:
            def __init__(self):
                self.headers = {'Content-Type': 'image/jpeg'}
                self.url = 'https://example.org/download'

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'\xff\xd8\xff\xe0fakejpeg'

            def close(self):
                return None

        class JPEGSession:
            def get(self, *args, **kwargs):
                return JPEGResponse()

        result = download_media(JPEGSession(), 'https://example.org/download', str(destination))

        assert result['status'] == 'downloaded'
        assert result['destination_path'].endswith('.jpg')
        assert os.path.exists(result['destination_path'])
        assert not os.path.exists(destination)

    def test_direct_download_uses_umask_adjusted_mode(self, tmp_path):
        destination = tmp_path / 'out' / 'image.bin'

        class JPEGResponse:
            headers = {'Content-Type': 'image/jpeg'}
            url = 'https://example.org/download'

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'\xff\xd8\xff\xe0fakejpeg'

            def close(self):
                return None

        class JPEGSession:
            def get(self, *args, **kwargs):
                return JPEGResponse()

        previous_umask = os.umask(0o027)
        try:
            result = download_media(
                JPEGSession(),
                'https://example.org/download',
                str(destination),
            )
        finally:
            os.umask(previous_umask)

        if os.name != 'nt':
            assert stat.S_IMODE(os.stat(result['destination_path']).st_mode) == 0o640

    def test_download_media_reuses_existing_cache_variant_extension(self, tmp_path):
        destination = tmp_path / 'out' / 'image.bin'
        cache_path = tmp_path / 'cache' / 'image.bin'
        cache_variant = tmp_path / 'cache' / 'image.jpg'
        cache_variant.parent.mkdir(parents=True)
        cache_variant.write_bytes(b'\xff\xd8\xff\xe0jpeg-data')

        class FailingSession:
            def get(self, *args, **kwargs):
                raise AssertionError('network should not be used when cache variant exists')

        result = download_media(FailingSession(), 'https://example.org/download', str(destination), cache_path=str(cache_path))

        assert result['status'] == 'cached'
        assert result['destination_path'].endswith('.jpg')
        assert (tmp_path / 'out' / 'image.jpg').read_bytes().startswith(b'\xff\xd8\xff')

    def test_download_media_rejects_html_returned_from_image_url(self, tmp_path):
        destination = tmp_path / 'out' / 'error.jpg'

        class HTMLResponse:
            headers = {'Content-Type': 'text/html'}
            url = 'https://example.org/error.jpg'

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                yield b'<html><body>rate limited</body></html>'

            def close(self):
                return None

        class HTMLSession:
            def get(self, *args, **kwargs):
                return HTMLResponse()

        with pytest.raises(MediaDownloadError, match='non-image Content-Type'):
            download_media(HTMLSession(), 'https://example.org/error.jpg', str(destination))

        assert not destination.exists()
        assert list(destination.parent.glob('*.tmp')) == []

    def test_download_media_enforces_content_length_limit(self, tmp_path):
        destination = tmp_path / 'out' / 'large.jpg'

        class LargeResponse:
            headers = {'Content-Type': 'image/jpeg', 'Content-Length': '1000'}

            def raise_for_status(self):
                return None

            def iter_content(self, chunk_size=65536):
                raise AssertionError('oversized response should be rejected before streaming')

            def close(self):
                return None

        class LargeSession:
            def get(self, *args, **kwargs):
                return LargeResponse()

        with pytest.raises(MediaDownloadError, match='exceeding --max-download-bytes'):
            download_media(
                LargeSession(),
                'https://example.org/large.jpg',
                str(destination),
                max_download_bytes=100,
            )

        assert not destination.exists()
        assert list(destination.parent.glob('*.tmp')) == []


class TestImagePostprocessing:
    def test_processed_raster_uses_umask_and_preserves_existing_mode(
        self,
        tmp_path,
    ):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'image.png'
        destination = tmp_path / 'image.jpg'
        Image.new('RGB', (40, 20), 'black').save(source)

        previous_umask = os.umask(0o027)
        try:
            result = postprocess_media_file(
                str(source),
                make_image_args(output_format='jpg', max_edge=16),
            )
        finally:
            os.umask(previous_umask)
        assert result == str(destination)
        if os.name != 'nt':
            assert stat.S_IMODE(destination.stat().st_mode) == 0o640

        Image.new('RGB', (40, 20), 'black').save(source)
        destination.chmod(0o664)
        postprocess_media_file(
            str(source),
            make_image_args(output_format='jpg', max_edge=12),
        )
        if os.name != 'nt':
            assert stat.S_IMODE(destination.stat().st_mode) == 0o664

    def test_postprocess_media_file_resizes_and_pads_raster(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'image.png'
        Image.new('RGB', (40, 20), 'black').save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                max_edge=16,
                canvas='square',
                background='white',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (16, 16)
        assert result.endswith('.png')

    def test_postprocess_media_file_trims_white_border(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'trim.png'
        image = Image.new('RGB', (20, 20), 'white')
        for x in range(5, 15):
            for y in range(7, 13):
                image.putpixel((x, y), (0, 0, 0))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='white',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (10, 6)

    def test_postprocess_media_file_trim_shape_square_center_crops_trimmed_result(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'trim-square.png'
        image = Image.new('RGBA', (24, 24), (0, 0, 0, 0))
        for x in range(8, 14):
            for y in range(4, 20):
                image.putpixel((x, y), (0, 0, 0, 255))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='transparent',
                trim_shape='square',
                background='transparent',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (6, 6)
            assert processed.mode == 'RGBA'
            assert processed.getchannel('A').getextrema() == (255, 255)

    def test_postprocess_media_file_trims_semantic_largest_rgb_component(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'semantic-rgb.png'
        image = Image.new('RGB', (40, 24), 'white')
        for x in range(3, 9):
            for y in range(2, 20):
                image.putpixel((x, y), (0, 0, 0))
        for x in range(20, 35):
            for y in range(5, 17):
                image.putpixel((x, y), (0, 128, 0))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='semantic',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (15, 12)
            assert processed.getpixel((0, 0)) == (0, 128, 0)

    def test_postprocess_media_file_trims_semantic_largest_alpha_component(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'semantic-alpha.png'
        image = Image.new('RGBA', (30, 24), (0, 0, 0, 0))
        for x in range(2, 8):
            for y in range(3, 21):
                image.putpixel((x, y), (255, 0, 0, 255))
        for x in range(14, 26):
            for y in range(6, 16):
                image.putpixel((x, y), (0, 0, 255, 255))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='semantic',
                trim_shape='square',
                background='transparent',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (10, 10)
            assert processed.mode == 'RGBA'
            assert processed.getpixel((0, 0)) == (0, 0, 255, 255)

    def test_postprocess_media_file_trim_shape_square_without_trim_center_crops_full_image(self, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'full-square.png'
        image = Image.new('RGB', (20, 10))
        for x in range(20):
            for y in range(10):
                image.putpixel((x, y), (x, 0, 0))
        image.save(source)

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                trim='off',
                trim_shape='square',
                background='white',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (10, 10)
            assert processed.getpixel((0, 0)) == (5, 0, 0)
            assert processed.getpixel((9, 0)) == (14, 0, 0)

    def test_postprocess_media_file_rasterizes_svg_when_requested(self, monkeypatch, tmp_path):
        Image = pytest.importorskip('PIL.Image')
        source = tmp_path / 'shape.svg'
        source.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"></svg>')

        monkeypatch.setattr(
            'nwkit.image.rasterize_svg_to_image',
            lambda source_path, max_edge=None: Image.new('RGBA', (20, 10), (0, 0, 0, 255)),
        )

        result = postprocess_media_file(
            str(source),
            make_image_args(
                output_format='png',
                max_edge=12,
                canvas='square',
                background='white',
            ),
        )

        with Image.open(result) as processed:
            assert processed.size == (12, 12)
        assert result.endswith('.png')

    def test_postprocess_media_file_requires_cairosvg_for_svg_conversion(self, monkeypatch, tmp_path):
        source = tmp_path / 'shape.svg'
        source.write_text('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"></svg>')

        monkeypatch.setattr(
            'nwkit.image.load_cairosvg_module',
            lambda: (_ for _ in ()).throw(RuntimeError('SVG image post-processing requires the optional CairoSVG dependency. Install optional image-processing dependencies with: pip install "nwkit[image]"')),
        )

        with pytest.raises(RuntimeError, match='CairoSVG dependency'):
            postprocess_media_file(
                str(source),
                make_image_args(
                    output_format='png',
                ),
            )

    def test_cairosvg_native_library_error_has_install_guidance(self, monkeypatch):
        real_import = builtins.__import__

        def import_with_missing_cairo(name, *args, **kwargs):
            if name == 'cairosvg':
                raise OSError('no library called cairo was found')
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, '__import__', import_with_missing_cairo)

        with pytest.raises(RuntimeError, match='native Cairo library') as error:
            image_module.load_cairosvg_module()
        assert 'pip install -e ".[image]"' in str(error.value)

    def test_image_main_skips_optional_processing_deps_when_not_requested(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'

        phylopic_candidates = {
            'Apis mellifera': [{
                'provider': 'phylopic',
                'provider_record_id': 'phy-plain',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'PhyloPic Artist',
                'source_page_url': 'https://api.phylopic.org/images/phy-plain',
                'media_url': 'https://images.phylopic.org/apis.svg',
                'width': 1000,
                'height': 900,
                'asset_type': 'silhouette',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'phylopic': DummyProvider(phylopic_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            with open(destination_path, 'wb') as handle:
                handle.write(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')
            return {'status': 'downloaded', 'destination_path': destination_path}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)
        monkeypatch.setattr(
            'nwkit.image.load_pillow_modules',
            lambda: (_ for _ in ()).throw(AssertionError('Pillow should not be loaded without post-processing options')),
        )

        image_main(
            make_image_args(
                infile=str(tree_path),
                out_dir=str(out_dir),
                source='phylopic',
            )
        )

        manifest_rows = read_tsv(out_dir / 'manifest.tsv')
        assert manifest_rows[0]['local_path'].endswith('.svg')

    def test_image_main_requires_pillow_when_postprocessing_is_requested(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'

        inaturalist_candidates = {
            'Apis mellifera': [{
                'provider': 'inaturalist',
                'provider_record_id': 'inat-plain',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Photographer',
                'source_page_url': 'https://www.inaturalist.org/observations/1',
                'media_url': 'https://static.inaturalist.org/apis.jpg',
                'width': 1000,
                'height': 900,
                'asset_type': 'photo',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'inaturalist': DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            write_valid_test_media(destination_path)
            return {'status': 'downloaded', 'destination_path': destination_path}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)
        monkeypatch.setattr(
            'nwkit.image.load_pillow_modules',
            lambda: (_ for _ in ()).throw(RuntimeError('Image post-processing requires the optional Pillow dependency. Install optional image-processing dependencies with: pip install "nwkit[image]"')),
        )

        with pytest.raises(RuntimeError, match='Pillow dependency'):
            image_main(
                make_image_args(
                    infile=str(tree_path),
                    out_dir=str(out_dir),
                    source='inaturalist',
                    output_format='png',
                )
            )

    def test_image_main_requires_pillow_when_trim_shape_square_is_requested(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'

        phylopic_candidates = {
            'Apis mellifera': [{
                'provider': 'phylopic',
                'provider_record_id': 'phy-plain',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'PhyloPic Artist',
                'source_page_url': 'https://api.phylopic.org/images/phy-plain',
                'media_url': 'https://images.phylopic.org/apis.svg',
                'width': 1000,
                'height': 900,
                'asset_type': 'silhouette',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'phylopic': DummyProvider(phylopic_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            with open(destination_path, 'wb') as handle:
                handle.write(b'<svg xmlns="http://www.w3.org/2000/svg"></svg>')
            return {'status': 'downloaded', 'destination_path': destination_path}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)
        monkeypatch.setattr('nwkit.image.rasterize_svg_to_image', lambda source_path, max_edge=None: object())
        monkeypatch.setattr(
            'nwkit.image.load_pillow_modules',
            lambda: (_ for _ in ()).throw(RuntimeError('Image post-processing requires the optional Pillow dependency. Install optional image-processing dependencies with: pip install "nwkit[image]"')),
        )

        with pytest.raises(RuntimeError, match='Pillow dependency'):
            image_main(
                make_image_args(
                    infile=str(tree_path),
                    out_dir=str(out_dir),
                    source='phylopic',
                    trim_shape='square',
                )
            )

    def test_image_main_requires_pillow_when_trim_semantic_is_requested(self, monkeypatch, tmp_path):
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Apis_mellifera_A);')
        out_dir = tmp_path / 'out'

        inaturalist_candidates = {
            'Apis mellifera': [{
                'provider': 'inaturalist',
                'provider_record_id': 'inat-semantic',
                'matched_name': 'Apis mellifera',
                'matched_rank': 'species',
                'license_code': 'cc-by',
                'license_url': 'https://creativecommons.org/licenses/by/4.0/',
                'attribution': 'Photographer',
                'source_page_url': 'https://www.inaturalist.org/observations/1',
                'media_url': 'https://static.inaturalist.org/apis.jpg',
                'width': 1000,
                'height': 900,
                'asset_type': 'photo',
            }],
        }

        def fake_build_providers(args, sources, session=None):
            providers = {
                'inaturalist': DummyProvider(inaturalist_candidates),
            }
            return DummySession(), None, providers

        def fake_download_media(session, media_url, destination_path, cache_path=None, max_download_bytes=None, provider=None, **kwargs):
            write_valid_test_media(destination_path)
            return {'status': 'downloaded', 'destination_path': destination_path}

        monkeypatch.setattr('nwkit.image.build_providers', fake_build_providers)
        monkeypatch.setattr('nwkit.image.download_media', fake_download_media)
        monkeypatch.setattr(
            'nwkit.image.load_pillow_modules',
            lambda: (_ for _ in ()).throw(RuntimeError('Image post-processing requires the optional Pillow dependency. Install optional image-processing dependencies with: pip install "nwkit[image]"')),
        )

        with pytest.raises(RuntimeError, match='Pillow dependency'):
            image_main(
                make_image_args(
                    infile=str(tree_path),
                    out_dir=str(out_dir),
                    source='inaturalist',
                    trim='semantic',
                )
            )


class TestImageSecurityLimits:
    @pytest.mark.parametrize(
        'allowed_hosts',
        [None, ('attacker.example',)],
        ids=['arbitrary-host', 'allowlisted-host'],
    )
    def test_real_session_pins_the_single_validated_dns_result(
        self,
        monkeypatch,
        allowed_hosts,
    ):
        resolver_calls = []
        pinned_calls = []
        public_address = '93.184.216.34'

        def alternating_resolution(hostname):
            resolver_calls.append(hostname)
            if len(resolver_calls) == 1:
                return {public_address}
            return {'127.0.0.1'}

        class Response:
            status_code = 200
            headers = {}
            url = 'https://attacker.example/image.png'

            def close(self):
                return None

        def fake_pinned_request(session, method, url, address, **kwargs):
            pinned_calls.append((method, url, address))
            return Response()

        monkeypatch.setattr(
            image_module,
            'resolve_hostname_addresses',
            alternating_resolution,
        )
        monkeypatch.setattr(
            image_module,
            '_request_pinned_address',
            fake_pinned_request,
        )
        session = requests.Session()
        try:
            response = image_module.safe_external_request(
                session,
                'get',
                'https://attacker.example/image.png',
                allowed_hosts=allowed_hosts,
            )
        finally:
            session.close()

        assert response.status_code == 200
        assert resolver_calls == ['attacker.example']
        assert pinned_calls == [
            ('get', 'https://attacker.example/image.png', public_address)
        ]

    def test_pinned_adapter_connects_to_ip_with_original_tls_hostname(self):
        adapter = image_module._PinnedHTTPSAdapter(
            address='93.184.216.34',
            tls_hostname='example.com',
        )
        request = requests.Request('GET', 'https://example.com/image.png').prepare()
        try:
            pool = adapter.get_connection_with_tls_context(
                request,
                verify=True,
            )
            assert pool.host == '93.184.216.34'
            assert pool.conn_kw['server_hostname'] == 'example.com'
            assert pool.assert_hostname == 'example.com'
        finally:
            adapter.close()

    def test_pinned_request_rejects_proxy_configuration(self):
        session = requests.Session()
        session.trust_env = False
        session.proxies = {'https': 'http://proxy.example:8080'}
        try:
            with pytest.raises(MediaDownloadError, match='Proxy'):
                image_module._request_pinned_address(
                    session,
                    'get',
                    'https://example.com/image.png',
                    '93.184.216.34',
                )
        finally:
            session.close()

    def test_cross_origin_redirect_drops_session_and_request_credentials(
        self,
        monkeypatch,
    ):
        calls = list()

        class Response:
            def __init__(self, status_code, url, location=None):
                self.status_code = status_code
                self.url = url
                self.headers = {} if location is None else {'Location': location}

            def close(self):
                return None

        responses = iter([
            Response(
                302,
                'https://origin.example/start',
                'https://other.example/next',
            ),
            Response(200, 'https://other.example/next'),
        ])

        def fake_get(pinned_session, url, **kwargs):
            effective_headers = {
                str(key).lower(): value
                for key, value in pinned_session.headers.items()
            }
            effective_headers.update({
                str(key).lower(): value
                for key, value in (kwargs.get('headers') or {}).items()
            })
            calls.append({
                'url': url,
                'auth': pinned_session.auth,
                'cert': pinned_session.cert,
                'headers': effective_headers,
                'cookies': pinned_session.cookies.get_dict(),
                'params': dict(pinned_session.params),
                'request_auth': kwargs.get('auth'),
                'request_cert': kwargs.get('cert'),
                'request_cookies': kwargs.get('cookies'),
                'request_params': kwargs.get('params'),
            })
            return next(responses)

        monkeypatch.setattr(requests.Session, 'get', fake_get)
        monkeypatch.setattr(
            image_module,
            '_resolve_public_addresses',
            lambda hostname, url: ('93.184.216.34',),
        )
        session = requests.Session()
        session.trust_env = False
        session.auth = ('session-user', 'session-secret')
        session.cert = '/tmp/session-client.pem'
        session.params = {'session-api-key': 'secret'}
        session.headers.update({
            'Authorization': 'Bearer session-secret',
            'Cookie': 'session-cookie=secret',
            'Proxy-Authorization': 'Basic proxy-secret',
            'X-Api-Key': 'session-api-secret',
        })
        session.cookies.set('ambient-cookie', 'secret')
        try:
            response = image_module.safe_external_request(
                session,
                'get',
                'https://origin.example/start',
                headers={
                    'Authorization': 'Bearer request-secret',
                    'Cookie': 'request-cookie=secret',
                    'Proxy-Authorization': 'Basic request-proxy-secret',
                    'X-Api-Key': 'request-api-secret',
                    'Accept-Language': 'ja',
                },
                auth=('request-user', 'request-secret'),
                cert='/tmp/request-client.pem',
                cookies={'request-cookie': 'secret'},
                params={'request-api-key': 'secret'},
            )
            response.close()
        finally:
            session.close()

        assert calls[0]['auth'] == ('session-user', 'session-secret')
        assert calls[0]['cert'] == '/tmp/session-client.pem'
        assert calls[0]['params'] == {'session-api-key': 'secret'}
        assert calls[0]['request_auth'] == ('request-user', 'request-secret')
        assert calls[0]['request_cert'] == '/tmp/request-client.pem'
        assert calls[0]['request_cookies'] == {'request-cookie': 'secret'}
        assert calls[0]['request_params'] == {'request-api-key': 'secret'}
        assert calls[0]['cookies'] == {'ambient-cookie': 'secret'}
        assert calls[1]['auth'] is None
        assert calls[1]['cert'] is None
        assert calls[1]['params'] == {}
        assert calls[1]['request_auth'] is None
        assert calls[1]['request_cert'] is None
        assert calls[1]['request_cookies'] is None
        assert calls[1]['request_params'] is None
        assert calls[1]['cookies'] == {}
        assert calls[1]['headers']['accept-language'] == 'ja'
        assert 'authorization' not in calls[1]['headers']
        assert 'cookie' not in calls[1]['headers']
        assert 'proxy-authorization' not in calls[1]['headers']
        assert 'x-api-key' not in calls[1]['headers']

    def test_pinned_redirect_merges_response_cookies_into_caller_session(
        self,
        monkeypatch,
    ):
        cookies_seen = list()

        class Response:
            def __init__(self, status_code, url, location=None):
                self.status_code = status_code
                self.url = url
                self.headers = {} if location is None else {'Location': location}

            def close(self):
                return None

        def fake_get(pinned_session, url, **kwargs):
            cookies_seen.append(pinned_session.cookies.get_dict(
                domain='cookie.example',
                path='/',
            ))
            if len(cookies_seen) == 1:
                pinned_session.cookies.set(
                    'route',
                    'blue',
                    domain='cookie.example',
                    path='/',
                )
                return Response(302, url, '/next')
            return Response(200, url)

        monkeypatch.setattr(requests.Session, 'get', fake_get)
        monkeypatch.setattr(
            image_module,
            '_resolve_public_addresses',
            lambda hostname, url: ('93.184.216.34',),
        )
        session = requests.Session()
        session.trust_env = False
        try:
            response = image_module.safe_external_request(
                session,
                'get',
                'https://cookie.example/start',
            )
            response.close()
            assert session.cookies.get(
                'route',
                domain='cookie.example',
                path='/',
            ) == 'blue'
        finally:
            session.close()

        assert cookies_seen == [{}, {'route': 'blue'}]

    def test_custom_session_without_redirect_control_is_not_retried(self):
        calls = list()

        class LegacySession:
            def get(self, url, **kwargs):
                calls.append(dict(kwargs))
                if 'allow_redirects' in kwargs:
                    raise TypeError(
                        "get() got an unexpected keyword argument 'allow_redirects'"
                    )
                raise AssertionError('unsafe fallback request was attempted')

        with pytest.raises(MediaDownloadError, match='allow_redirects=False'):
            image_module.safe_external_request(
                LegacySession(),
                'get',
                'https://public.example/start',
            )

        assert len(calls) == 1
        assert calls[0]['allow_redirects'] is False

    def test_mixed_public_and_private_dns_answers_are_rejected(self, monkeypatch):
        monkeypatch.setattr(
            image_module,
            'resolve_hostname_addresses',
            lambda hostname: {'93.184.216.34', '127.0.0.1'},
        )
        monkeypatch.setattr(
            image_module,
            '_request_pinned_address',
            lambda *args, **kwargs: pytest.fail('request must not be sent'),
        )
        session = requests.Session()
        try:
            with pytest.raises(MediaDownloadError, match='non-public address'):
                image_module.safe_external_request(
                    session,
                    'get',
                    'https://mixed.example/image.png',
                )
        finally:
            session.close()

    def test_rejects_private_redirect_before_following_it(self):
        calls = list()

        class RedirectResponse:
            status_code = 302
            headers = {'Location': 'https://127.0.0.1/private'}
            url = 'https://example.org/start'

            def close(self):
                return None

        class RedirectSession:
            def get(self, url, **kwargs):
                calls.append(url)
                return RedirectResponse()

        with pytest.raises(MediaDownloadError, match='non-public address'):
            image_module.safe_external_request(
                RedirectSession(),
                'get',
                'https://example.org/start',
            )
        assert calls == ['https://example.org/start']

    def test_download_rejects_provider_host_outside_allowlist(self, tmp_path):
        class UnusedSession:
            def get(self, *args, **kwargs):
                raise AssertionError('request should not be sent')

        with pytest.raises(MediaDownloadError, match='not allowed'):
            download_media(
                UnusedSession(),
                'https://example.org/silhouette.svg',
                str(tmp_path / 'image.svg'),
                provider='phylopic',
            )

    def test_provider_json_body_is_bounded_while_streaming(self):
        class LargeResponse:
            headers = {}

            def iter_content(self, chunk_size):
                yield b'{"value":"'
                yield b'x' * 32
                yield b'"}'

        with pytest.raises(MediaDownloadError, match='exceeded'):
            image_module.bounded_response_json(LargeResponse(), max_bytes=16)

    def test_cached_candidate_with_insecure_url_is_discarded(self, tmp_path):
        cache_path = tmp_path / 'query.json'
        image_module.write_cached_provider_candidates(
            str(cache_path),
            [{
                'provider': 'openverse',
                'media_url': 'http://example.org/image.jpg',
            }],
            fetch_limit=10,
        )
        assert image_module.load_cached_provider_candidates(str(cache_path)) is None

    @pytest.mark.parametrize(
        'loader',
        [
            image_module.load_cached_provider_candidates,
            image_module.load_cached_bioicons_catalog,
        ],
    )
    def test_non_object_json_cache_is_treated_as_a_cache_miss(
        self,
        loader,
        tmp_path,
    ):
        cache_path = tmp_path / 'cache.json'
        cache_path.write_text('[]')

        assert loader(str(cache_path)) is None

    def test_default_svg_output_rejects_external_content(self, tmp_path):
        source = tmp_path / 'unsafe.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<image href="https://example.org/tracker.png"/>'
            '</svg>'
        )
        with pytest.raises(MediaDownloadError, match='forbidden|external'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_svg_output_allows_safe_inline_and_style_tag_css(self, tmp_path):
        source = tmp_path / 'safe-css.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<defs>'
            '<filter id="shadow"><feOffset dx="1" dy="1"/></filter>'
            '<style>.safe { fill:red; stroke:#000; filter:url(#shadow); }</style>'
            '</defs>'
            '<rect class="safe" width="20" height="10" '
            'style="fill:blue;stroke:red"/>'
            '</svg>'
        )

        result = postprocess_media_file(str(source), make_image_args())

        assert result == str(source)

    @pytest.mark.parametrize(
        'unsafe_css',
        [
            '<rect style="fill:url(https://attacker.example/pattern.svg)"/>',
            '<rect style="fill:url(images/pattern.svg)"/>',
            '<style>@import "https://attacker.example/theme.css";</style>',
            '<rect style="width:expression(alert(1))"/>',
            '<rect style="fill:javascript:alert(1)"/>',
        ],
        ids=[
            'absolute-url',
            'relative-url',
            'import',
            'expression',
            'javascript',
        ],
    )
    def test_default_svg_output_rejects_unsafe_css_references(
        self,
        tmp_path,
        unsafe_css,
    ):
        source = tmp_path / 'unsafe-css.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '{}'
            '</svg>'.format(unsafe_css)
        )

        with pytest.raises(MediaDownloadError, match='external|executable'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_svg_output_rejects_css_escape_obfuscation(self, tmp_path):
        source = tmp_path / 'escaped-url.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<rect width="20" height="10" '
            'style="fill:u\\72l(https://attacker.example/p.svg#x)"/>'
            '</svg>'
        )

        with pytest.raises(MediaDownloadError, match='external|executable'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_raster_output_rejects_truncated_png(self, tmp_path):
        source = tmp_path / 'truncated.png'
        source.write_bytes(b'\x89PNG\r\n\x1a\n\x00\x00\x00\x0d')

        with pytest.raises(MediaDownloadError, match='truncated|corrupt|unsupported'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_raster_output_rejects_truncated_jpeg(self, tmp_path):
        source = tmp_path / 'truncated.jpg'
        Image, _, _ = image_module.load_pillow_modules()
        Image.new('RGB', (64, 64), 'red').save(source, format='JPEG')
        source.write_bytes(source.read_bytes()[:-2])

        with pytest.raises(MediaDownloadError, match='truncated|corrupt|unsupported'):
            postprocess_media_file(str(source), make_image_args())

    def test_svg_rejects_utf16_before_entity_parsing(self, tmp_path):
        source = tmp_path / 'utf16-entity.svg'
        source.write_bytes((
            '<?xml version="1.0" encoding="UTF-16"?>'
            '<!DOCTYPE svg [<!ENTITY local "expanded">]>'
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<text>&local;</text>'
            '</svg>'
        ).encode('utf-16'))

        with pytest.raises(MediaDownloadError, match='UTF-16|encoding'):
            postprocess_media_file(str(source), make_image_args())

    def test_default_raster_output_rejects_oversized_dimensions(self, tmp_path):
        source = tmp_path / 'oversized.png'
        Image, _, _ = image_module.load_pillow_modules()
        Image.new('RGBA', (1, 1), (0, 0, 0, 255)).save(source)
        payload = bytearray(source.read_bytes())
        payload[16:24] = struct.pack('>II', 100_000, 100_000)
        payload[29:33] = struct.pack(
            '>I',
            zlib.crc32(bytes(payload[12:29])) & 0xffffffff,
        )
        source.write_bytes(payload)

        with pytest.raises(MediaDownloadError) as error:
            postprocess_media_file(str(source), make_image_args())
        assert 'pixel' in str(error.value).lower() or 'bomb' in str(error.value).lower()

    def test_default_svg_output_strips_external_doctype(self, tmp_path):
        source = tmp_path / 'doctype.svg'
        source.write_text(
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" '
            '"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"/>'
        )
        source.chmod(0o664)

        result = postprocess_media_file(str(source), make_image_args())

        assert result == str(source)
        assert '<!DOCTYPE' not in source.read_text()
        if os.name != 'nt':
            assert stat.S_IMODE(source.stat().st_mode) == 0o664

    def test_svg_dimensions_are_read_from_sanitized_xml_root(
        self,
        monkeypatch,
        tmp_path,
    ):
        source = tmp_path / 'sanitized-dimensions.svg'
        source.write_text(
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE svg SYSTEM "https://example.org/external.dtd">\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10"/>'
        )
        monkeypatch.setattr(
            image_module.ElementTree,
            'parse',
            lambda *args, **kwargs: pytest.fail(
                'validated SVG dimensions must not reparse the original file'
            ),
        )

        result = postprocess_media_file(str(source), make_image_args())

        assert result == str(source)
        assert '<!DOCTYPE' not in source.read_text()

    def test_rasterizing_local_svg_does_not_modify_input(self, tmp_path):
        pytest.importorskip('cairosvg')
        source = tmp_path / 'local-input.svg'
        source.write_text(
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" '
            '"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">\n'
            '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="10">'
            '<rect width="20" height="10" fill="#000"/>'
            '</svg>'
        )
        before = source.read_bytes()

        rasterized = image_module.rasterize_svg_to_image(str(source))

        assert rasterized.size == (20, 10)
        assert source.read_bytes() == before

    def test_svg_is_scaled_before_rasterization(self, tmp_path):
        pytest.importorskip('cairosvg')
        pytest.importorskip('PIL.Image')
        source = tmp_path / 'huge.svg'
        source.write_text(
            '<svg xmlns="http://www.w3.org/2000/svg" width="100000" height="50000">'
            '<rect width="100000" height="50000" fill="#000"/>'
            '</svg>'
        )
        rasterized = image_module.rasterize_svg_to_image(str(source), max_edge=100)
        assert rasterized.size == (100, 50)

    def test_ncbi_extraction_rejects_oversized_member(self, monkeypatch, tmp_path):
        archive_path = tmp_path / 'new_taxdump.tar.gz'
        destination_path = tmp_path / 'images.dmp'
        payload = b'12345'
        with tarfile.open(archive_path, 'w:gz') as archive:
            info = tarfile.TarInfo('images.dmp')
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
        monkeypatch.setattr(image_module, 'NCBI_IMAGES_TABLE_MAX_BYTES', 4)

        with pytest.raises(MediaDownloadError, match='outside the supported range'):
            image_module._extract_ncbi_images_table(
                str(archive_path),
                str(destination_path),
            )
        assert not destination_path.exists()

    def test_ncbi_cache_replaces_corrupt_images_table(self, monkeypatch, tmp_path):
        args = make_image_args(
            download_dir=str(tmp_path / 'downloads'),
            out_dir=str(tmp_path / 'out'),
        )
        cache_dir = tmp_path / 'downloads' / 'ncbi-taxonomy-images'
        cache_dir.mkdir(parents=True)
        images_path = cache_dir / 'images.dmp'
        images_path.write_text('corrupt cache')
        valid_row = (
            b'64365\t|\timage:Cyanophora paradoxa\t|\t'
            b'https://example.org/image.jpg\t|\tCC BY 4.0\t|\tAuthor\t|\t'
            b'Publisher\t|\t\t|\t2762\t|\n'
        )

        def fake_download(session, url, destination_path, **kwargs):
            with tarfile.open(destination_path, 'w:gz') as archive:
                info = tarfile.TarInfo('images.dmp')
                info.size = len(valid_row)
                archive.addfile(info, io.BytesIO(valid_row))

        monkeypatch.setattr(image_module, '_fetch_ncbi_archive_md5', lambda session: None)
        monkeypatch.setattr(image_module, '_download_to_path', fake_download)

        result = image_module.ensure_ncbi_images_table(args, DummySession())

        assert result == str(images_path)
        assert images_path.read_bytes() == valid_row

    def test_ncbi_database_validation_handles_uri_characters_and_empty_tables(
        self,
        tmp_path,
    ):
        cache_dir = tmp_path / ('cache #%' if os.name == 'nt' else 'cache?#%')
        cache_dir.mkdir()
        database_path = cache_dir / 'images.sqlite3'
        connection = sqlite3.connect(database_path)
        connection.execute(
            'CREATE TABLE images ('
            'taxid INTEGER, record_id TEXT, image_url TEXT)'
        )
        connection.commit()
        connection.close()

        assert not image_module._ncbi_images_database_is_valid(
            str(database_path)
        )

        connection = sqlite3.connect(database_path)
        connection.execute(
            'INSERT INTO images VALUES (?, ?, ?)',
            (1, 'record', 'https://example.org/image.jpg'),
        )
        connection.commit()
        connection.close()

        assert image_module._ncbi_images_database_is_valid(str(database_path))

    def test_ncbi_database_reloads_table_after_late_corrupt_row(
        self,
        monkeypatch,
        tmp_path,
    ):
        args = make_image_args(
            download_dir=str(tmp_path / 'downloads'),
            out_dir=str(tmp_path / 'out'),
        )
        cache_dir = tmp_path / 'downloads' / 'ncbi-taxonomy-images'
        images_path = cache_dir / 'images.dmp'
        valid_row = (
            '64365\t|\timage:Cyanophora paradoxa\t|\t'
            'https://example.org/image.jpg\t|\tCC BY 4.0\t|\tAuthor\t|\t'
            'Publisher\t|\t\t|\t2762\t|\n'
        )
        calls = {'count': 0}

        def fake_ensure_table(args, session):
            calls['count'] += 1
            cache_dir.mkdir(parents=True, exist_ok=True)
            payload = (
                valid_row + 'corrupt later row\n'
                if calls['count'] == 1
                else valid_row
            )
            images_path.write_text(payload)
            return str(images_path)

        monkeypatch.setattr(
            image_module,
            'ensure_ncbi_images_table',
            fake_ensure_table,
        )

        database_path = image_module.ensure_ncbi_images_database(
            args,
            DummySession(),
        )

        assert calls['count'] == 2
        assert image_module._ncbi_images_database_is_valid(database_path)

    def test_image_outputs_must_be_distinct(self, tmp_path):
        out_dir = tmp_path / 'out'
        with pytest.raises(ValueError, match='Output paths must be distinct'):
            image_main(make_image_args(
                out_dir=str(out_dir),
                manifest_out=str(out_dir / 'ATTRIBUTION.md'),
                source='phylopic',
            ))
