"""Unit tests for geeutils.py."""
import pytest
import unittest.mock as mock

import ee
import hydrafloods as hf


class TestUtils:
    def test_export_image_toAsset(self):
        """Test that the toAsset export gets called by default."""
        # mock ee export functionality
        ee.batch.Export.image.toAsset = mock.MagicMock()

        # call export function
        hf.geeutils.export_image(ee.Image(), ee.Geometry.Point(0, 0), "assid", "desc")
        # assert that the ee.export call within the function was called
        assert ee.batch.Export.image.toAsset.call_count == 1

    def test_export_image_toDrive(self):
        """Test that the toDrive export gets called when specified."""
        # mock ee export functionality
        ee.batch.Export.image.toDrive = mock.MagicMock()

        # call export function
        hf.geeutils.export_image(
            ee.Image(), ee.Geometry.Point(0, 0), "assid", "desc", export_type="toDrive"
        )
        # assert that the ee.export call within the function was called
        assert ee.batch.Export.image.toDrive.call_count == 1

    def test_export_image_export_type_invalid_string(self):
        """Test ValueError for export_type."""
        # mock ee export functionality
        ee.batch.Export.image.toAsset = mock.MagicMock()

        # call export function
        with pytest.raises(ValueError):
            hf.geeutils.export_image(
                ee.Image(),
                ee.Geometry.Point(0, 0),
                "assid",
                "desc",
                export_type="invalid string",
            )

    def test_export_image_export_type_invalid_type(self):
        """Test TypeError for export_type."""
        # mock ee export functionality
        ee.batch.Export.image.toAsset = mock.MagicMock()

        # call export function
        with pytest.raises(TypeError):
            hf.geeutils.export_image(
                ee.Image(), ee.Geometry.Point(0, 0), "assid", "desc", export_type=123
            )

    def test_export_image_folder_invalid_type(self):
        """Test TypeError for folder."""
        # mock ee export functionality
        ee.batch.Export.image.toAsset = mock.MagicMock()

        # call export function
        with pytest.raises(TypeError):
            hf.geeutils.export_image(
                ee.Image(), ee.Geometry.Point(0, 0), "assid", "desc", folder=123
            )
